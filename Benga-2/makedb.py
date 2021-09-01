import os
import logging
import argparse
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict, Counter

import pandas as pd
from matplotlib import pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline

from utils import run_cmd
from profiling import profiling, filter_blast_result
plt.style.use('ggplot')

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
console = logging.StreamHandler()
console.setLevel(logging.INFO)


def parallel(fun, args_ls, threads):
    with Pool(threads) as p:
        try:
            for args in args_ls:
                p.apply_async(fun, args=args)
            p.close()
            p.join()
        except KeyboardInterrupt:
            p.terminate()


def remove_non_cds_feature(in_gff, out_gff):
    with open(in_gff) as in_f, open(out_gff, 'w') as out_f:
        for line in in_f:
            if line != '##FASTA\n':
                if line.startswith('##'):
                    out_f.write(line)
                else:
                    if line.split()[2] == 'CDS':
                        out_f.write(line)
            else:
                out_f.write(line)
                out_f.write(in_f.read())


def annotate_genome(genome, out_path, prefix, prodigaltf):
    cmd = f"prokka --prefix {prefix} --cpus 1 --outdir {out_path} {genome}"
    if prodigaltf:
        cmd += f" --prodigaltf {prodigaltf}"
    run_cmd(cmd)


def plot_genome_coverage(data, figure_path, core=95):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(data, bins=range(102), histtype='step', lw=2)
    ax.set_title("Genome coverage distribution", fontsize=22)
    ax.set_xlabel("Percentage of genomes covered by loci (%)", fontsize=18)
    ax.set_ylabel("Number of locus", fontsize=18)
    ax.axvline(core, ls='--', color='b')
    fig.savefig(figure_path, dpi=300)


def workflow(file_list, output_path, threads, prodigaltf):
    root_path = Path(output_path)
    annotation_path = root_path/'Annotated'

    logging.info("Annotated")
    args_ls = []
    for file in file_list:
        file_path = Path(file)
        args = (file_path, annotation_path/file_path.stem, file_path.stem, prodigaltf)
        args_ls.append(args)
    parallel(annotate_genome, args_ls, threads)

    gff_path = root_path/'GFF'
    gff_path.mkdir()
    for dirpath in annotation_path.iterdir():
        src_gff = dirpath/(dirpath.name + '.gff')
        dst_gff = gff_path/(dirpath.name + '.gff')
        remove_non_cds_feature(src_gff, dst_gff)

    logging.info("Pan genome analysis")
    roary_path = root_path/'roary'
    run_cmd(f"roary -p {threads} -i 95 -s -f {roary_path} {gff_path/'*.gff'}")

    clustered_proteins = roary_path / "clustered_proteins"
    allele_to_locus_map = dict()
    f = open(clustered_proteins)
    for line in f:
        locus_tag, allele_tags = line.split(": ")
        allele_to_locus_map.update({allele_tag: locus_tag for allele_tag in allele_tags.split()})
    f.close()

    allele_frequency_count = defaultdict(Counter)
    for ffn_file in annotation_path.glob('**/*.ffn'):
        for record in SeqIO.parse(ffn_file, 'fasta'):
            if record.id in allele_to_locus_map:
                locus_tag = allele_to_locus_map[record.id]
                allele_frequency_count[locus_tag].update([record.seq])

    logging.info("Define locus reference sequence")
    pan_genome_nucl_sequence = root_path/'pan_genome.fna'
    pan_genome_prot_sequence = root_path/'pan_genome.faa'
    with open(pan_genome_nucl_sequence, 'w') as nf, open(pan_genome_prot_sequence, 'w') as pf:
        for locus_tag, count in allele_frequency_count.items():
            representative = count.most_common(1)[0][0]
            n_record = SeqRecord(representative, id=locus_tag, description='')
            p_record = SeqRecord(representative.translate(table=11), id=locus_tag, description='')
            SeqIO.write(n_record, nf, 'fasta')
            SeqIO.write(p_record, pf, 'fasta')

    logging.info("Profiling")
    profile_path = root_path/'Profile'
    profile_path.mkdir()

    args_ls = []
    for dirpath in annotation_path.iterdir():
        genome = dirpath / (dirpath.name + '.fna')
        outfile = profile_path / (dirpath.name + '.tsv')
        args = (genome, pan_genome_prot_sequence, outfile, prodigaltf, 1)
        args_ls.append(args)
    parallel(profiling, args_ls, threads)

    logging.info("Recalculate loci occurrence")
    loci_freq_count = Counter()
    num_isolate = 0
    for filepath in profile_path.iterdir():
        profile = pd.read_csv(filepath, sep='\t', usecols=['locus_id'])
        loci_freq_count.update(profile['locus_id'])
        num_isolate += 1
    occurrence = {locus: number/num_isolate*100 for locus, number in loci_freq_count.items()}

    gene_presence_absence = pd.read_csv(
        roary_path/'gene_presence_absence.csv',
        usecols=['Gene', 'Annotation', 'No. isolates'],
        index_col=0
    )
    new_gene_presence_absence = gene_presence_absence.filter(occurrence.keys(), axis=0)
    new_gene_presence_absence['No. isolates'] = new_gene_presence_absence.index.map(loci_freq_count)
    new_gene_presence_absence['occurrence'] = new_gene_presence_absence.index.map(occurrence)

    logging.info('Drop duplicate loci')
    cline = NcbimakeblastdbCommandline(
        input_file=pan_genome_prot_sequence, dbtype='prot',
    )
    cline(stdout=False)

    blast_out = root_path/'self_align'
    cline = NcbiblastpCommandline(
        query=pan_genome_prot_sequence,
        db=pan_genome_prot_sequence,
        evalue=1e-6, outfmt='6 qseqid sseqid pident length qlen slen', num_threads=threads,
        out=blast_out
    )
    cline(stdout=False)
    run_cmd(f"rm {os.path.join(root_path, 'pan_genome.*')}")

    drop_locus = set()
    for qseqid, sseqid in filter_blast_result(blast_out):
        if occurrence[qseqid] > occurrence[sseqid]:
            drop_locus.add(sseqid)
        else:
            drop_locus.add(qseqid)
    blast_out.unlink()

    pan_genome_results = new_gene_presence_absence.drop(index=drop_locus)
    pan_genome_results['Nucleotide'] = pan_genome_results.index.map({record.id: str(record.seq) for record in records})
    pan_genome_results['Peptide'] = pan_genome_results.index.map(
        {record.id: str(record.seq.translate(table=11)) for record in records})
    pan_genome_results.to_csv(root_path/'pan_genome_info.txt', sep='\t')

    logging.info('Plot genome coverage')
    plot_genome_coverage(
        pan_genome_results['occurrence'],
        root_path/'genome_coverage.png',
    )
    plot_genome_coverage(
        pan_genome_results[pan_genome_results['occurrence'] > 5]['occurrence'],
        root_path/'genome_coverage_5_prec.png',
    )
    logging.info('Done')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        nargs='+',
                        required=True,
                        help="Path of query genomes.")
    parser.add_argument("-o", "--output_path",
                        required=True,
                        help="Path of output directory.")
    parser.add_argument("--prodigaltf",
                        default='',
                        help="Path of prodigal training file. default:''")
    parser.add_argument("-t", "--threads",
                        type=int,
                        default=2,
                        help="Number of threads. default: 2")
    args = parser.parse_args()

    workflow(args.input, args.output_path, args.threads, args.prodigaltf)


if __name__ == '__main__':
    main()

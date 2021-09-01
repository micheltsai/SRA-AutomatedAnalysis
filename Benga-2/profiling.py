import os
import hashlib
import argparse
from tempfile import TemporaryDirectory

import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline, NcbimakeblastdbCommandline

from utils import run_cmd


def prodigal_cli(bin='prodigal', **kwargs):
    """
    Build prodigal command line.
    """
    command = [bin]
    parameters = {
        'a': '-a',
        'c': '-c',
        'd': '-d',
        'f': '-f',
        'g': '-g',
        'i': '-i',
        'm': '-m',
        'n': '-n',
        'o': '-o',
        'p': '-p',
        'q': '-q',
        's': '-s',
        't': '-t',
    }
    for key, value in kwargs.items():
        parameter = parameters.get(key)
        if parameter:
            if isinstance(value, bool):
                command += [parameter]
            elif value:
                command += [parameter, str(value)]
            else:
                pass
        else:
            raise ValueError('Parameter %s is invalided variable' % key)
    return ' '.join(command)


def make_blast_database(input_file, out, dbtype='prot'):
    cline = NcbimakeblastdbCommandline(input_file=input_file, dbtype=dbtype, out=out)
    cline()


def encode_sequence_id(record):
    """
    Convert nucleotide sequence to hash codes by hash function 'sha256'.
    """
    return hashlib.sha256(str(record.seq).encode("ascii")).hexdigest()


def predict_open_reading_frame(infile, outfile, prodigaltf=''):
    """
    Predict open reading frame and export nucleotide sequence
    """
    cmd = prodigal_cli(i=infile, d=outfile, g=11, c=True, m=True, q=True, t=prodigaltf)
    run_cmd(cmd)


def blastp(query_records, db, out, threads=2):
    query = (record.translate(table=11, id=True).format('fasta') for record in query_records)
    cline = NcbiblastpCommandline(
        out=out,
        db=db,
        evalue=1e-6,
        num_threads=threads,
        outfmt='6 qseqid sseqid pident length qlen slen',
    )
    cline(stdin=''.join(query))


def filter_blast_result(blast_out, identity=95, min_cov=.75, max_cov=1.25):
    """
    Filter identity >= 95, ratio of query sequence length to subject sequence length between 0.75 and 1.25
    """
    handle = open(blast_out)
    for line in handle:
        qseqid, sseqid, pident, length, qlen, slen = line.strip().split()
        pident, length, qlen, slen = float(pident), float(length), float(qlen), float(slen)
        if qseqid != sseqid and pident >= identity and min_cov <= qlen/slen < max_cov and min_cov <= qlen/length < max_cov:
            yield sseqid, qseqid
    handle.close()


def profiling(infile, scheme, outfile, prodigaltf, threads=2, allele_nucleotide_output=False):
    with TemporaryDirectory(dir="/tmp") as tmpdir:
        nucl_frames = os.path.join(tmpdir, 'frames.fna')
        predict_open_reading_frame(infile, nucl_frames, prodigaltf=prodigaltf)

        blast_db = os.path.join(tmpdir, 'scheme')
        make_blast_database(scheme, blast_db)

        nucl_records = []
        for record in SeqIO.parse(nucl_frames, 'fasta'):
            record.id = encode_sequence_id(record)
            nucl_records.append(record)

        blast_output = os.path.join(tmpdir, 'blast.out')
        blastp(nucl_records, blast_db, blast_output, threads)

        df = pd.DataFrame(filter_blast_result(blast_output), columns=['locus_id', 'allele_id'])
        df = df.sort_values('allele_id', kind='mergesort').drop_duplicates('locus_id')
        df.to_csv(outfile, sep='\t', index=False)

    if allele_nucleotide_output:
        hits = df['allele_id'].to_dict()
        f = open(allele_nucleotide_output, 'w')
        for record in nucl_records:
            if record.id in hits:
                record.description = hits[record.id]
                f.write(record.format('fasta'))
        f.close()


def main():
    parser = argparse.ArgumentParser('Benga')
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Path of query genome.")
    parser.add_argument("-o", "--output",
                        required=True,
                        help="Path of output file.")
    parser.add_argument("--scheme",
                        required=True,
                        help="Core-genome MLST scheme, it should is peptide sequence.")
    parser.add_argument("--prodigaltf",
                        default='',
                        help="Path of prodigal training file. default: None")
    parser.add_argument("-t", "--threads",
                        type=int,
                        default=2,
                        help="Number of threads. default: 2")
    parser.add_argument("--allele-nucleotide-output",
                        default="",
                        help="Output nucleotide FASTA file of allele nucleotide sequences. default: ''")
    args = parser.parse_args()

    profiling(
        infile=args.input,
        scheme=args.scheme,
        outfile=args.output_path,
        prodigaltf=args.prodigaltf,
        threads=args.threads,
        allele_nucleotide_output=args.allele_nucleotide_output
    )


if __name__ == '__main__':
    main()

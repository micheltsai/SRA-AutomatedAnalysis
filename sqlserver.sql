-- Exported from QuickDBD: https://www.quickdatabasediagrams.com/
-- Link to schema: https://app.quickdatabasediagrams.com/#/d/dBamRU
-- NOTE! If you have used non-SQL datatypes in your design, you will have to change these here.

-- Modify this code to update the DB schema diagram.
-- To reset the sample schema, replace everything with
-- two dots ('..' - without quotes).

SET XACT_ABORT ON

BEGIN TRANSACTION QUICKDBD

CREATE TABLE [MLST] (
    [mId] int  NOT NULL ,
    [SRAID] varchar(200)  NOT NULL ,
    [Profile] varchar(300)  NOT NULL ,
    [Organism] varchar(300)  NOT NULL ,
    [SequenceType] int  NOT NULL ,
    CONSTRAINT [PK_MLST] PRIMARY KEY CLUSTERED (
        [SRAID] ASC
    ),
    CONSTRAINT [UK_MLST_mId] UNIQUE (
        [mId]
    )
)

CREATE TABLE [Sistr] (
    [sId] int  NOT NULL ,
    [SRAID] varchar(200)  NOT NULL ,
    [Serovar] varchar(300)  NOT NULL ,
    CONSTRAINT [PK_Sistr] PRIMARY KEY CLUSTERED (
        [SRAID] ASC
    ),
    CONSTRAINT [UK_Sistr_sId] UNIQUE (
        [sId]
    )
)

CREATE TABLE [Amrfinder] (
    [aId] int  NOT NULL ,
    [SRAID] varchar(200)  NOT NULL ,
    [GeneSymbol] varchar(300)  NOT NULL ,
    [ElementSubtype] varchaar(300)  NOT NULL ,
    [Method] varchar(200)  NOT NULL ,
    CONSTRAINT [PK_Amrfinder] PRIMARY KEY CLUSTERED (
        [SRAID] ASC
    )
)

-- Table documentation comment 1 (try the PDF/RTF export)
-- Table documentation comment 2
CREATE TABLE [SRA] (
    [Id] int  NOT NULL ,
    [Genome] varchar(200)  NOT NULL ,
    [Seq] varchar(300)  NOT NULL ,
    CONSTRAINT [PK_SRA] PRIMARY KEY CLUSTERED (
        [Id] ASC
    ),
    CONSTRAINT [UK_SRA_Genome] UNIQUE (
        [Genome]
    )
)

CREATE TABLE [Plasmidfinder] (
    [PlasId] int  NOT NULL ,
    [SRAID] varchar(200)  NOT NULL ,
    [pDB] varchar(300)  NOT NULL ,
    [Plasmid] varchar(200)  NOT NULL ,
    CONSTRAINT [PK_Plasmidfinder] PRIMARY KEY CLUSTERED (
        [SRAID] ASC
    ),
    CONSTRAINT [UK_Plasmidfinder_PlasId] UNIQUE (
        [PlasId]
    )
)

CREATE TABLE [PlasDB] (
    [pDBid] int  NOT NULL ,
    [DBNmae] varchar(300)  NOT NULL ,
    CONSTRAINT [PK_PlasDB] PRIMARY KEY CLUSTERED (
        [pDBid] ASC
    )
)

CREATE TABLE [mlstOrganism] (
    [oid] int  NOT NULL ,
    [oName] varchar(300)  NOT NULL ,
    CONSTRAINT [PK_mlstOrganism] PRIMARY KEY CLUSTERED (
        [oid] ASC
    )
)

ALTER TABLE [MLST] WITH CHECK ADD CONSTRAINT [FK_MLST_SRAID] FOREIGN KEY([SRAID])
REFERENCES [SRA] ([Genome])

ALTER TABLE [MLST] CHECK CONSTRAINT [FK_MLST_SRAID]

ALTER TABLE [Sistr] WITH CHECK ADD CONSTRAINT [FK_Sistr_SRAID] FOREIGN KEY([SRAID])
REFERENCES [SRA] ([Genome])

ALTER TABLE [Sistr] CHECK CONSTRAINT [FK_Sistr_SRAID]

ALTER TABLE [Amrfinder] WITH CHECK ADD CONSTRAINT [FK_Amrfinder_SRAID] FOREIGN KEY([SRAID])
REFERENCES [SRA] ([Genome])

ALTER TABLE [Amrfinder] CHECK CONSTRAINT [FK_Amrfinder_SRAID]

ALTER TABLE [Plasmidfinder] WITH CHECK ADD CONSTRAINT [FK_Plasmidfinder_SRAID] FOREIGN KEY([SRAID])
REFERENCES [SRA] ([Genome])

ALTER TABLE [Plasmidfinder] CHECK CONSTRAINT [FK_Plasmidfinder_SRAID]

ALTER TABLE [PlasDB] WITH CHECK ADD CONSTRAINT [FK_PlasDB_DBNmae] FOREIGN KEY([DBNmae])
REFERENCES [Plasmidfinder] ([pDB])

ALTER TABLE [PlasDB] CHECK CONSTRAINT [FK_PlasDB_DBNmae]

ALTER TABLE [mlstOrganism] WITH CHECK ADD CONSTRAINT [FK_mlstOrganism_oName] FOREIGN KEY([oName])
REFERENCES [MLST] ([Organism])

ALTER TABLE [mlstOrganism] CHECK CONSTRAINT [FK_mlstOrganism_oName]

COMMIT TRANSACTION QUICKDBD
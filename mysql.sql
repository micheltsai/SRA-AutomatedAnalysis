-- Exported from QuickDBD: https://www.quickdatabasediagrams.com/
-- Link to schema: https://app.quickdatabasediagrams.com/#/d/dBamRU
-- NOTE! If you have used non-SQL datatypes in your design, you will have to change these here.

-- Modify this code to update the DB schema diagram.
-- To reset the sample schema, replace everything with
-- two dots ('..' - without quotes).

CREATE TABLE `MLST` (
    `mId` int  NOT NULL ,
    `SRAID` varchar(200)  NOT NULL ,
    `Profile` varchar(300)  NOT NULL ,
    `Organism` varchar(300)  NOT NULL ,
    `SequenceType` int  NOT NULL ,
    PRIMARY KEY (
        `SRAID`
    ),
    CONSTRAINT `uc_MLST_mId` UNIQUE (
        `mId`
    )
);

CREATE TABLE `Sistr` (
    `sId` int  NOT NULL ,
    `SRAID` varchar(200)  NOT NULL ,
    `Serovar` varchar(300)  NOT NULL ,
    PRIMARY KEY (
        `SRAID`
    ),
    CONSTRAINT `uc_Sistr_sId` UNIQUE (
        `sId`
    )
);

CREATE TABLE `Amrfinder` (
    `aId` int  NOT NULL ,
    `SRAID` varchar(200)  NOT NULL ,
    `GeneSymbol` varchar(300)  NOT NULL ,
    `ElementSubtype` varchaar(300)  NOT NULL ,
    `Method` varchar(200)  NOT NULL ,
    PRIMARY KEY (
        `SRAID`
    )
);

-- Table documentation comment 1 (try the PDF/RTF export)
-- Table documentation comment 2
CREATE TABLE `SRA` (
    `Id` int  NOT NULL ,
    `Genome` varchar(200)  NOT NULL ,
    `Seq` varchar(300)  NOT NULL ,
    PRIMARY KEY (
        `Id`
    ),
    CONSTRAINT `uc_SRA_Genome` UNIQUE (
        `Genome`
    )
);

CREATE TABLE `Plasmidfinder` (
    `PlasId` int  NOT NULL ,
    `SRAID` varchar(200)  NOT NULL ,
    `pDB` varchar(300)  NOT NULL ,
    `Plasmid` varchar(200)  NOT NULL ,
    PRIMARY KEY (
        `SRAID`
    ),
    CONSTRAINT `uc_Plasmidfinder_PlasId` UNIQUE (
        `PlasId`
    )
);

CREATE TABLE `PlasDB` (
    `pDBid` int  NOT NULL ,
    `DBNmae` varchar(300)  NOT NULL ,
    PRIMARY KEY (
        `pDBid`
    )
);

CREATE TABLE `mlstOrganism` (
    `oid` int  NOT NULL ,
    `oName` varchar(300)  NOT NULL ,
    PRIMARY KEY (
        `oid`
    )
);

ALTER TABLE `MLST` ADD CONSTRAINT `fk_MLST_SRAID` FOREIGN KEY(`SRAID`)
REFERENCES `SRA` (`Genome`);

ALTER TABLE `Sistr` ADD CONSTRAINT `fk_Sistr_SRAID` FOREIGN KEY(`SRAID`)
REFERENCES `SRA` (`Genome`);

ALTER TABLE `Amrfinder` ADD CONSTRAINT `fk_Amrfinder_SRAID` FOREIGN KEY(`SRAID`)
REFERENCES `SRA` (`Genome`);

ALTER TABLE `Plasmidfinder` ADD CONSTRAINT `fk_Plasmidfinder_SRAID` FOREIGN KEY(`SRAID`)
REFERENCES `SRA` (`Genome`);

ALTER TABLE `PlasDB` ADD CONSTRAINT `fk_PlasDB_DBNmae` FOREIGN KEY(`DBNmae`)
REFERENCES `Plasmidfinder` (`pDB`);

ALTER TABLE `mlstOrganism` ADD CONSTRAINT `fk_mlstOrganism_oName` FOREIGN KEY(`oName`)
REFERENCES `MLST` (`Organism`);


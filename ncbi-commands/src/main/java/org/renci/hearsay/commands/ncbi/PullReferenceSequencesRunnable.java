package org.renci.hearsay.commands.ncbi;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.renci.gene2accession.G2AFilter;
import org.renci.gene2accession.G2AParser;
import org.renci.gene2accession.filter.G2AAndFilter;
import org.renci.gene2accession.filter.G2AAssemblyFilter;
import org.renci.gene2accession.filter.G2AGenomicNucleotideAccessionVersionPrefixFilter;
import org.renci.gene2accession.filter.G2ARNANucleotideAccessionVersionPrefixFilter;
import org.renci.gene2accession.filter.G2ATaxonIdFilter;
import org.renci.gene2accession.model.OrientationType;
import org.renci.gene2accession.model.Record;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.model.Gene;
import org.renci.hearsay.dao.model.GenomeReference;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.ReferenceSequenceType;
import org.renci.hearsay.dao.model.StrandType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullReferenceSequencesRunnable implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(PullReferenceSequencesRunnable.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullReferenceSequencesRunnable() {
        super();
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {
            // 380MB gzipped
            File genes2RefSeqFile = FTPUtil.ncbiDownload("/gene/DATA", "gene2refseq.gz");
            G2AParser gene2AccessionParser = G2AParser.getInstance(8);
            G2AAndFilter andFilter = new G2AAndFilter(Arrays.asList(new G2AFilter[] { new G2ATaxonIdFilter(9606),
                    new G2AAssemblyFilter("Reference GRCh38.p2 Primary Assembly"),
                    new G2AGenomicNucleotideAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NC_" })),
                    new G2ARNANucleotideAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })) }));
            List<Record> recordList = gene2AccessionParser.parse(andFilter, genes2RefSeqFile);
            for (Record record : recordList) {

                if (StringUtils.isEmpty(record.getRNANucleotideAccessionVersion())) {
                    continue;
                }

                ReferenceSequence referenceSequence = new ReferenceSequence();
                referenceSequence
                        .setStrandType(record.getOrientation().equals(OrientationType.MINUS) ? StrandType.MINUS
                                : StrandType.PLUS);

                String prefix = record.getRNANucleotideAccessionVersion().substring(0, 3);
                for (ReferenceSequenceType referenceSequenceType : ReferenceSequenceType.values()) {
                    if (referenceSequenceType.getPrefixes().contains(prefix)) {
                        referenceSequence.setType(referenceSequenceType);
                        break;
                    }
                }

                Location genomicLocation = new Location(record.getGenomicStartPosition(),
                        record.getGenomicEndPosition());
                genomicLocation.setId(hearsayDAOBean.getLocationDAO().save(genomicLocation));
                referenceSequence.setGenomicLocation(genomicLocation);

                referenceSequence.setId(hearsayDAOBean.getReferenceSequenceDAO().save(referenceSequence));
                logger.info(referenceSequence.toString());

                // set Gene
                Gene exampleGene = new Gene();
                exampleGene.setSymbol(record.getSymbol());
                List<Gene> potentialGenes = hearsayDAOBean.getGeneDAO().findByExample(exampleGene);
                if (potentialGenes != null && !potentialGenes.isEmpty()) {
                    referenceSequence.setGene(potentialGenes.get(0));
                }

                // set GenomeReference
                GenomeReference exampleGenomeReference = new GenomeReference();
                exampleGenomeReference.setName(record.getAssembly().replace("Reference", "")
                        .replace("Primary Assembly", "").trim());
                List<GenomeReference> potentialGenomeReferences = hearsayDAOBean.getGenomeReferenceDAO().findByExample(
                        exampleGenomeReference);
                if (potentialGenomeReferences != null && !potentialGenomeReferences.isEmpty()) {
                    referenceSequence.setGenomeReference(potentialGenomeReferences.get(0));
                }

                // set nucleotide identifier
                String versionedRefSeqAccession = record.getRNANucleotideAccessionVersion();
                Identifier identifier = new Identifier("www.ncbi.nlm.nih.gov/nuccore", versionedRefSeqAccession);
                List<Identifier> possibleIdentifiers = hearsayDAOBean.getIdentifierDAO().findByExample(identifier);
                if (possibleIdentifiers != null && !possibleIdentifiers.isEmpty()) {
                    identifier = possibleIdentifiers.get(0);
                } else {
                    identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                }
                logger.debug(identifier.toString());
                referenceSequence.getIdentifiers().add(identifier);

                // set protein identifier
                String versionedProteinAccession = record.getProteinAccessionVersion();
                identifier = new Identifier("www.ncbi.nlm.nih.gov/protein", versionedProteinAccession);
                possibleIdentifiers = hearsayDAOBean.getIdentifierDAO().findByExample(identifier);
                if (possibleIdentifiers != null && !possibleIdentifiers.isEmpty()) {
                    identifier = possibleIdentifiers.get(0);
                } else {
                    identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                }
                logger.debug(identifier.toString());
                referenceSequence.getIdentifiers().add(identifier);

                // set genomic identifier
                String versionedGenomicAccession = record.getGenomicNucleotideAccessionVersion();
                identifier = new Identifier("www.ncbi.nlm.nih.gov/genome", versionedGenomicAccession);
                possibleIdentifiers = hearsayDAOBean.getIdentifierDAO().findByExample(identifier);
                if (possibleIdentifiers != null && !possibleIdentifiers.isEmpty()) {
                    identifier = possibleIdentifiers.get(0);
                } else {
                    identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                }
                logger.debug(identifier.toString());
                referenceSequence.getIdentifiers().add(identifier);

                hearsayDAOBean.getReferenceSequenceDAO().save(referenceSequence);

            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}

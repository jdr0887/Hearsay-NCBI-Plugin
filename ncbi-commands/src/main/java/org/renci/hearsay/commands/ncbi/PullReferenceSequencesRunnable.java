package org.renci.hearsay.commands.ncbi;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.collections4.CollectionUtils;
import org.renci.gene2accession.G2AFilter;
import org.renci.gene2accession.G2AParser;
import org.renci.gene2accession.filter.G2AAndFilter;
import org.renci.gene2accession.filter.G2AAssemblyFilter;
import org.renci.gene2accession.filter.G2AGenomicNucleotideAccessionVersionPrefixFilter;
import org.renci.gene2accession.filter.G2AProteinAccessionVersionPrefixFilter;
import org.renci.gene2accession.filter.G2ARNANucleotideAccessionVersionPrefixFilter;
import org.renci.gene2accession.filter.G2ATaxonIdFilter;
import org.renci.gene2accession.model.OrientationType;
import org.renci.gene2accession.model.Record;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBeanService;
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

    private static final Logger logger = LoggerFactory.getLogger(PullReferenceSequencesRunnable.class);

    private HearsayDAOBeanService hearsayDAOBeanService;

    public PullReferenceSequencesRunnable(HearsayDAOBeanService hearsayDAOBeanService) {
        super();
        this.hearsayDAOBeanService = hearsayDAOBeanService;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {
            // 380MB gzipped
            File genes2RefSeqFile = FTPUtil.ncbiDownload("/gene/DATA", "gene2refseq.gz");
            G2AParser gene2AccessionParser = G2AParser.getInstance(8);
            List<G2AFilter> filters = Arrays.asList(new G2AFilter[] { new G2ATaxonIdFilter(9606),
                    // new G2AAssemblyFilter("Reference.*(Primary Assembly|ALT_REF_LOCI.*)"),
                    new G2AAssemblyFilter("Reference.*Primary Assembly"),
                    new G2AProteinAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NP_" })),
                    new G2AGenomicNucleotideAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NC_" })),
                    new G2ARNANucleotideAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })) });
            G2AAndFilter andFilter = new G2AAndFilter(filters);
            List<Record> recordList = gene2AccessionParser.parse(andFilter, genes2RefSeqFile);
            ExecutorService es = Executors.newFixedThreadPool(8);
            for (Record record : recordList) {
                es.submit(new PersistReferenceSequencesRunnable(record));
            }
            es.shutdown();
            es.awaitTermination(1L, TimeUnit.HOURS);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    class PersistReferenceSequencesRunnable implements Runnable {

        private final Record record;

        public PersistReferenceSequencesRunnable(Record record) {
            super();
            this.record = record;
        }

        @Override
        public void run() {
            logger.debug("ENTERING run()");

            try {
                ReferenceSequence referenceSequence = new ReferenceSequence();
                referenceSequence.setStrandType(record.getOrientation().equals(OrientationType.MINUS) ? StrandType.MINUS : StrandType.PLUS);

                String prefix = record.getRNANucleotideAccessionVersion().substring(0, 3);
                for (ReferenceSequenceType referenceSequenceType : ReferenceSequenceType.values()) {
                    if (referenceSequenceType.getPrefixes().contains(prefix)) {
                        referenceSequence.setType(referenceSequenceType);
                        break;
                    }
                }

                Location genomicLocation = new Location(record.getGenomicStartPosition(), record.getGenomicEndPosition());
                genomicLocation.setId(hearsayDAOBeanService.getLocationDAO().save(genomicLocation));
                referenceSequence.setGenomicLocation(genomicLocation);

                referenceSequence.setId(hearsayDAOBeanService.getReferenceSequenceDAO().save(referenceSequence));

                // set Gene
                Gene exampleGene = new Gene();
                exampleGene.setSymbol(record.getSymbol());
                List<Gene> potentialGenes = hearsayDAOBeanService.getGeneDAO().findByExample(exampleGene);
                if (CollectionUtils.isNotEmpty(potentialGenes)) {
                    referenceSequence.setGene(potentialGenes.get(0));
                }

                // set GenomeReference
                GenomeReference exampleGenomeReference = new GenomeReference();
                exampleGenomeReference.setName(record.getAssembly().replace("Reference", "").replace("Primary Assembly", "").trim());
                List<GenomeReference> potentialGenomeReferences = hearsayDAOBeanService.getGenomeReferenceDAO()
                        .findByExample(exampleGenomeReference);
                if (CollectionUtils.isNotEmpty(potentialGenomeReferences)) {
                    referenceSequence.setGenomeReference(potentialGenomeReferences.get(0));
                }

                // set nucleotide identifier
                String versionedRefSeqAccession = record.getRNANucleotideAccessionVersion();
                Identifier identifier = new Identifier("www.ncbi.nlm.nih.gov/nuccore", versionedRefSeqAccession);
                List<Identifier> possibleIdentifiers = hearsayDAOBeanService.getIdentifierDAO().findByExample(identifier);
                if (CollectionUtils.isNotEmpty(possibleIdentifiers)) {
                    identifier = possibleIdentifiers.get(0);
                } else {
                    identifier.setId(hearsayDAOBeanService.getIdentifierDAO().save(identifier));
                }
                referenceSequence.getIdentifiers().add(identifier);

                // set protein identifier
                String versionedProteinAccession = record.getProteinAccessionVersion();
                identifier = new Identifier("www.ncbi.nlm.nih.gov/protein", versionedProteinAccession);
                possibleIdentifiers = hearsayDAOBeanService.getIdentifierDAO().findByExample(identifier);
                if (CollectionUtils.isNotEmpty(possibleIdentifiers)) {
                    identifier = possibleIdentifiers.get(0);
                } else {
                    identifier.setId(hearsayDAOBeanService.getIdentifierDAO().save(identifier));
                }
                referenceSequence.getIdentifiers().add(identifier);

                // set genomic identifier
                String versionedGenomicAccession = record.getGenomicNucleotideAccessionVersion();
                identifier = new Identifier("www.ncbi.nlm.nih.gov/genome", versionedGenomicAccession);
                possibleIdentifiers = hearsayDAOBeanService.getIdentifierDAO().findByExample(identifier);
                if (CollectionUtils.isNotEmpty(possibleIdentifiers)) {
                    identifier = possibleIdentifiers.get(0);
                } else {
                    identifier.setId(hearsayDAOBeanService.getIdentifierDAO().save(identifier));
                }
                referenceSequence.getIdentifiers().add(identifier);

                hearsayDAOBeanService.getReferenceSequenceDAO().save(referenceSequence);
                logger.info("refSeqAccession = {}, proteinAccession = {}, genomicAccession = {}", versionedRefSeqAccession,
                        versionedProteinAccession, versionedGenomicAccession);
            } catch (Exception e) {
                logger.error(e.getMessage(), e);
                e.printStackTrace();
            }

        }

    }

}

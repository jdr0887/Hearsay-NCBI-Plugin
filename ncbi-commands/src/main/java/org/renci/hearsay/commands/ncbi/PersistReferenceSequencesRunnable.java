package org.renci.hearsay.commands.ncbi;

import java.util.List;

import org.renci.gene2accession.model.OrientationType;
import org.renci.gene2accession.model.Record;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.HearsayDAOException;
import org.renci.hearsay.dao.model.Gene;
import org.renci.hearsay.dao.model.GenomeReference;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.ReferenceSequenceType;
import org.renci.hearsay.dao.model.StrandType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PersistReferenceSequencesRunnable implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(PersistReferenceSequencesRunnable.class);

    private final Record record;

    private final HearsayDAOBean hearsayDAOBean;

    public PersistReferenceSequencesRunnable(Record record, HearsayDAOBean hearsayDAOBean) {
        super();
        this.record = record;
        this.hearsayDAOBean = hearsayDAOBean;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {
            ReferenceSequence referenceSequence = new ReferenceSequence();
            referenceSequence.setStrandType(record.getOrientation().equals(OrientationType.MINUS) ? StrandType.MINUS
                    : StrandType.PLUS);

            String prefix = record.getRNANucleotideAccessionVersion().substring(0, 3);
            for (ReferenceSequenceType referenceSequenceType : ReferenceSequenceType.values()) {
                if (referenceSequenceType.getPrefixes().contains(prefix)) {
                    referenceSequence.setType(referenceSequenceType);
                    break;
                }
            }

            Location genomicLocation = new Location(record.getGenomicStartPosition(), record.getGenomicEndPosition());
            genomicLocation.setId(hearsayDAOBean.getLocationDAO().save(genomicLocation));
            referenceSequence.setGenomicLocation(genomicLocation);

            referenceSequence.setId(hearsayDAOBean.getReferenceSequenceDAO().save(referenceSequence));

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
            referenceSequence.getIdentifiers().add(identifier);

            hearsayDAOBean.getReferenceSequenceDAO().save(referenceSequence);
            logger.info("refSeqAccession = {}, proteinAccession = {}, genomicAccession = {}", versionedRefSeqAccession,
                    versionedProteinAccession, versionedGenomicAccession);
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            e.printStackTrace();
        }

    }

}

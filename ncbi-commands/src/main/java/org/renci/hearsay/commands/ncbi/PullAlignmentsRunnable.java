package org.renci.hearsay.commands.ncbi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.collections.CollectionUtils;
import org.renci.gbff.GBFFFilter;
import org.renci.gbff.GBFFManager;
import org.renci.gbff.filter.GBFFAndFilter;
import org.renci.gbff.filter.GBFFFeatureSourceOrganismNameFilter;
import org.renci.gbff.filter.GBFFFeatureTypeNameFilter;
import org.renci.gbff.filter.GBFFSequenceAccessionPrefixFilter;
import org.renci.gbff.filter.GBFFSourceOrganismNameFilter;
import org.renci.gbff.model.Feature;
import org.renci.gbff.model.Sequence;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullAlignmentsRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PullAlignmentsRunnable.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullAlignmentsRunnable() {
        super();
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {

            LinkedList<String> alignmentsLines = new LinkedList<String>();

            File alignmentsFile = FTPUtil.ncbiDownload("/refseq/H_sapiens/alignments",
                    "GCF_000001405.28_knownrefseq_alignments.gff3");
            try (FileInputStream fis = new FileInputStream(alignmentsFile);
                    InputStreamReader isr = new InputStreamReader(fis);
                    BufferedReader br = new BufferedReader(isr)) {
                String line;
                while ((line = br.readLine()) != null) {
                    alignmentsLines.add(line);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            alignmentsFile.delete();

            // this will take a while
            GBFFManager gbffMgr = GBFFManager.getInstance(1, true);

            List<GBFFFilter> filters = Arrays.asList(new GBFFFilter[] {
                    new GBFFSequenceAccessionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })),
                    new GBFFSourceOrganismNameFilter("Homo sapiens"),
                    new GBFFFeatureSourceOrganismNameFilter("Homo sapiens"), new GBFFFeatureTypeNameFilter("CDS"),
                    new GBFFFeatureTypeNameFilter("source") });

            GBFFAndFilter gbffFilter = new GBFFAndFilter(filters);

            List<File> fileList = FTPUtil.ncbiDownloadBySuffix("/refseq/H_sapiens/mRNA_Prot", "rna.gbff.gz");

            ExecutorService es = Executors.newFixedThreadPool(8);
            for (File f : fileList) {
                List<Sequence> sequenceList = gbffMgr.deserialize(gbffFilter, f);
                logger.debug("sequenceList.size(): {}", sequenceList.size());
                if (CollectionUtils.isNotEmpty(sequenceList)) {

                    for (Sequence sequence : sequenceList) {

                        String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ") ? sequence
                                .getVersion().substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();

                        List<Identifier> rnaNucleotideAccessionIdentifierList = hearsayDAOBean
                                .getIdentifierDAO()
                                .findByExample(new Identifier("www.ncbi.nlm.nih.gov/nuccore", refSeqVersionedAccession));

                        String proteinAccession = null;
                        Feature firstCDSFeature = null;
                        for (Feature feature : sequence.getFeatures()) {
                            if (!"CDS".equals(feature.getType())) {
                                continue;
                            }
                            firstCDSFeature = feature;
                            break;
                        }
                        proteinAccession = firstCDSFeature.getQualifiers().getProperty("protein_id").replace("\"", "");

                        List<Identifier> proteinAccessionIdentifierList = hearsayDAOBean.getIdentifierDAO()
                                .findByExample(new Identifier("www.ncbi.nlm.nih.gov/protein", proteinAccession));

                        List<Long> identifierIdList = new ArrayList<Long>();
                        for (Identifier identifier : rnaNucleotideAccessionIdentifierList) {
                            identifierIdList.add(identifier.getId());
                        }
                        for (Identifier identifier : proteinAccessionIdentifierList) {
                            identifierIdList.add(identifier.getId());
                        }

                        List<ReferenceSequence> potentialRefSeqs = hearsayDAOBean.getReferenceSequenceDAO()
                                .findByIdentifiers(identifierIdList);

                        if (CollectionUtils.isEmpty(potentialRefSeqs)) {
                            logger.warn(
                                    "Could not find ReferenceSequence: refSeqVersionedAccession = {}, proteinAccession = {}",
                                    refSeqVersionedAccession, proteinAccession);
                            continue;
                        }

                        logger.info("Using ReferenceSequence: refSeqVersionedAccession = {}, proteinAccession = {}",
                                refSeqVersionedAccession, proteinAccession);

                        // persistAlignmentsExecutorService.submit(new PersistAlignmentsFromUCSCRunnable(hearsayDAOBean,
                        // sequence));
                        es.submit(new PersistAlignmentsFromNCBIRunnable(hearsayDAOBean, sequence, potentialRefSeqs,
                                firstCDSFeature, alignmentsLines));
                        es.submit(new PersistFeaturesRunnable(hearsayDAOBean, sequence, potentialRefSeqs));
                    }

                }
                // f.delete();
            }
            es.shutdown();
            es.awaitTermination(4L, TimeUnit.HOURS);

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }

    }

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}

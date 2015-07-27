package org.renci.hearsay.commands.ncbi;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.renci.gbff.GBFFFilter;
import org.renci.gbff.GBFFManager;
import org.renci.gbff.filter.GBFFAndFilter;
import org.renci.gbff.filter.GBFFFeatureSourceOrganismNameFilter;
import org.renci.gbff.filter.GBFFFeatureTypeNameFilter;
import org.renci.gbff.filter.GBFFSequenceAccessionPrefixFilter;
import org.renci.gbff.filter.GBFFSourceOrganismNameFilter;
import org.renci.gbff.model.Sequence;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBean;
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

            List<String> accessionPrefixList = Arrays.asList(new String[] { "NM_", "NR_" });

            GBFFManager gbffMgr = GBFFManager.getInstance(1, true);

            List<GBFFFilter> filters = Arrays.asList(new GBFFFilter[] {
                    new GBFFSequenceAccessionPrefixFilter(accessionPrefixList),
                    new GBFFSourceOrganismNameFilter("Homo sapiens"),
                    new GBFFFeatureSourceOrganismNameFilter("Homo sapiens"), new GBFFFeatureTypeNameFilter("CDS"),
                    new GBFFFeatureTypeNameFilter("source") });

            GBFFAndFilter gbffFilter = new GBFFAndFilter(filters);
            List<File> vertebrateMammalianFileList = FTPUtil.ncbiDownloadBySuffix(
                    "/refseq/release/vertebrate_mammalian", "rna.gbff.gz");

            for (File f : vertebrateMammalianFileList) {
                List<Sequence> sequenceList = gbffMgr.deserialize(gbffFilter, f);
                logger.info("sequenceList.size(): {}", sequenceList.size());
                if (sequenceList != null && !sequenceList.isEmpty()) {

                    ExecutorService persistAlignmentsExecutorService = Executors.newFixedThreadPool(4);
                    for (Sequence sequence : sequenceList) {
                        // persistAlignmentsExecutorService.submit(new PersistAlignmentsFromUCSCRunnable(hearsayDAOBean,
                        // sequence));
                        persistAlignmentsExecutorService.submit(new PersistAlignmentsFromNCBIRunnable(hearsayDAOBean,
                                sequence));
                    }
                    persistAlignmentsExecutorService.shutdown();
                    persistAlignmentsExecutorService.awaitTermination(1L, TimeUnit.HOURS);

                    ExecutorService persistFeaturesExecutorService = Executors.newFixedThreadPool(4);
                    for (Sequence sequence : sequenceList) {
                        persistFeaturesExecutorService.submit(new PersistFeaturesRunnable(hearsayDAOBean, sequence));
                    }
                    persistFeaturesExecutorService.shutdown();
                    persistFeaturesExecutorService.awaitTermination(1L, TimeUnit.HOURS);

                }
                // f.delete();
            }

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

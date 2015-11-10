package org.renci.hearsay.commands.ncbi;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.collections4.CollectionUtils;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.model.Chromosome;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PullRunnable.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullRunnable() {
        super();
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        // persist dictionary items
        long startPersistChromosomeTime = System.currentTimeMillis();
        try {
            List<String> chromosomeList = new ArrayList<String>();
            for (int i = 1; i < 22; i++) {
                chromosomeList.add(i + "");
            }
            chromosomeList.add("X");
            chromosomeList.add("Y");
            chromosomeList.add("MT");

            for (String chromosome : chromosomeList) {
                List<Chromosome> foundChromosomeList = hearsayDAOBean.getChromosomeDAO().findByName(chromosome);
                if (CollectionUtils.isEmpty(foundChromosomeList)) {
                    hearsayDAOBean.getChromosomeDAO().save(new Chromosome(chromosome));
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endPersistChromosomeTime = System.currentTimeMillis();

        // persist genes & genome references
        long startPersistGenesAndGenomeReferencesTime = System.currentTimeMillis();
        try {
            ExecutorService es = Executors.newFixedThreadPool(2);

            PullGenesRunnable pullGenesRunnable = new PullGenesRunnable();
            pullGenesRunnable.setHearsayDAOBean(hearsayDAOBean);
            es.submit(pullGenesRunnable);

            PullGenomeReferencesRunnable pullGenomeReferencesRunnable = new PullGenomeReferencesRunnable();
            pullGenomeReferencesRunnable.setHearsayDAOBean(hearsayDAOBean);
            es.submit(pullGenesRunnable);

            es.shutdown();
            es.awaitTermination(1L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long endPersistGenesAndGenomeReferencesTime = System.currentTimeMillis();

        long startPersistReferenceSequencesTime = System.currentTimeMillis();
        // persist reference sequences
        try {
            ExecutorService es = Executors.newSingleThreadExecutor();
            PullReferenceSequencesRunnable pullReferenceSequencesRunnable = new PullReferenceSequencesRunnable();
            pullReferenceSequencesRunnable.setHearsayDAOBean(hearsayDAOBean);
            es.submit(pullReferenceSequencesRunnable);
            es.shutdown();
            es.awaitTermination(2L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long endPersistReferenceSequencesTime = System.currentTimeMillis();

        // persist alignments
        long startPersistAlignmentsTime = System.currentTimeMillis();
        try {
            ExecutorService es = Executors.newSingleThreadExecutor();
            PullAlignmentsRunnable pullAlignmentsRunnable = new PullAlignmentsRunnable();
            pullAlignmentsRunnable.setHearsayDAOBean(hearsayDAOBean);
            es.submit(pullAlignmentsRunnable);
            es.shutdown();
            es.awaitTermination(6L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long endPersistAlignmentsTime = System.currentTimeMillis();

        // add alignment UTRs
        long startPersistAlignmentUTRsTime = System.currentTimeMillis();
        try {
            ExecutorService es = Executors.newSingleThreadExecutor();
            AddAlignmentUTRsRunnable addAlignmentUTRsRunnable = new AddAlignmentUTRsRunnable();
            addAlignmentUTRsRunnable.setHearsayDAOBean(hearsayDAOBean);
            es.submit(addAlignmentUTRsRunnable);
            es.shutdown();
            es.awaitTermination(2L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long endPersistAlignmentUTRsTime = System.currentTimeMillis();

        logger.info("duration to persist Chromosomes: {} seconds",
                (endPersistChromosomeTime - startPersistChromosomeTime) / 1000);

        logger.info("duration to persist Genes & GenomeReferences: {} seconds",
                (endPersistGenesAndGenomeReferencesTime - startPersistGenesAndGenomeReferencesTime) / 1000);

        logger.info("duration to persist ReferenceSequences: {} seconds",
                (endPersistReferenceSequencesTime - startPersistReferenceSequencesTime) / 1000);

        logger.info("duration to persist Alignments: {} seconds",
                (endPersistAlignmentsTime - startPersistAlignmentsTime) / 1000);

        logger.info("duration to persist Alignment UTRs: {} seconds",
                (endPersistAlignmentUTRsTime - startPersistAlignmentUTRsTime) / 1000);

    }

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}

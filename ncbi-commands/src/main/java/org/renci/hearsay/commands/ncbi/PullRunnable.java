package org.renci.hearsay.commands.ncbi;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.collections4.CollectionUtils;
import org.renci.hearsay.commands.AddAlignmentUTRsRunnable;
import org.renci.hearsay.dao.HearsayDAOBeanService;
import org.renci.hearsay.dao.model.Chromosome;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PullRunnable.class);

    private HearsayDAOBeanService hearsayDAOBeanService;

    public PullRunnable(HearsayDAOBeanService hearsayDAOBeanService) {
        super();
        this.hearsayDAOBeanService = hearsayDAOBeanService;
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
                List<Chromosome> foundChromosomes = hearsayDAOBeanService.getChromosomeDAO().findByName(chromosome);
                if (CollectionUtils.isEmpty(foundChromosomes)) {
                    hearsayDAOBeanService.getChromosomeDAO().save(new Chromosome(chromosome));
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endPersistChromosomeTime = System.currentTimeMillis();
        logger.info("duration to persist Chromosomes: {} seconds", (endPersistChromosomeTime - startPersistChromosomeTime) / 1000);

        // persist genes & genome references
        long startPersistGenesAndGenomeReferencesTime = System.currentTimeMillis();
        try {
            ExecutorService es = Executors.newFixedThreadPool(2);

            PullGenesRunnable pullGenesRunnable = new PullGenesRunnable(hearsayDAOBeanService);
            es.submit(pullGenesRunnable);

            PullGenomeReferencesRunnable pullGenomeReferencesRunnable = new PullGenomeReferencesRunnable(hearsayDAOBeanService);
            es.submit(pullGenomeReferencesRunnable);

            es.shutdown();
            es.awaitTermination(1L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long endPersistGenesAndGenomeReferencesTime = System.currentTimeMillis();
        logger.info("duration to persist Genes & GenomeReferences: {} seconds",
                (endPersistGenesAndGenomeReferencesTime - startPersistGenesAndGenomeReferencesTime) / 1000);

        // persist reference sequences
        long startPersistReferenceSequencesTime = System.currentTimeMillis();
        try {
            ExecutorService es = Executors.newSingleThreadExecutor();
            PullReferenceSequencesRunnable pullReferenceSequencesRunnable = new PullReferenceSequencesRunnable(hearsayDAOBeanService);
            es.submit(pullReferenceSequencesRunnable);
            es.shutdown();
            es.awaitTermination(2L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long endPersistReferenceSequencesTime = System.currentTimeMillis();
        logger.info("duration to persist ReferenceSequences: {} seconds",
                (endPersistReferenceSequencesTime - startPersistReferenceSequencesTime) / 1000);
        
        // persist alignments
        long startPersistAlignmentsTime = System.currentTimeMillis();
        try {
            ExecutorService es = Executors.newSingleThreadExecutor();
            PullAlignmentsRunnable pullAlignmentsRunnable = new PullAlignmentsRunnable(hearsayDAOBeanService);
            es.submit(pullAlignmentsRunnable);
            es.shutdown();
            es.awaitTermination(6L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long endPersistAlignmentsTime = System.currentTimeMillis();
        logger.info("duration to persist Alignments: {} seconds", (endPersistAlignmentsTime - startPersistAlignmentsTime) / 1000);

        // add alignment UTRs
        // long startPersistAlignmentUTRsTime = System.currentTimeMillis();
        // try {
        // ExecutorService es = Executors.newSingleThreadExecutor();
        // AddAlignmentUTRsRunnable addAlignmentUTRsRunnable = new AddAlignmentUTRsRunnable(hearsayDAOBeanService);
        // es.submit(addAlignmentUTRsRunnable);
        // es.shutdown();
        // es.awaitTermination(2L, TimeUnit.HOURS);
        // } catch (InterruptedException e) {
        // e.printStackTrace();
        // }
        // long endPersistAlignmentUTRsTime = System.currentTimeMillis();

        // logger.info("duration to persist Alignment UTRs: {} seconds", (endPersistAlignmentUTRsTime -
        // startPersistAlignmentUTRsTime) / 1000);

    }

}

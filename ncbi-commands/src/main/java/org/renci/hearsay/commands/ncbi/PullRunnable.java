package org.renci.hearsay.commands.ncbi;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.collections4.CollectionUtils;
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

        // persist reference sequences
        long startPersistReferenceSequencesTime = System.currentTimeMillis();
        try {
            PullReferenceSequencesRunnable pullReferenceSequencesRunnable = new PullReferenceSequencesRunnable(hearsayDAOBeanService);
            Executors.newSingleThreadExecutor().submit(pullReferenceSequencesRunnable).get();
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        long endPersistReferenceSequencesTime = System.currentTimeMillis();

        // persist alignments
        long startPersistAlignmentsTime = System.currentTimeMillis();
        try {
            PullAlignmentsRunnable pullAlignmentsRunnable = new PullAlignmentsRunnable(hearsayDAOBeanService);
            Executors.newSingleThreadExecutor().submit(pullAlignmentsRunnable).get();
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        long endPersistAlignmentsTime = System.currentTimeMillis();

        long startAddAlignmentUTRsTime = System.currentTimeMillis();
        try {
            AddAlignmentUTRsRunnable runnable = new AddAlignmentUTRsRunnable(hearsayDAOBeanService);
            Executors.newSingleThreadExecutor().submit(runnable).get();
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        long endAddAlignmentUTRsTime = System.currentTimeMillis();

        long startPersistFeaturesTime = System.currentTimeMillis();
        try {
            PullFeaturesRunnable pullFeaturesRunnable = new PullFeaturesRunnable(hearsayDAOBeanService);
            Executors.newSingleThreadExecutor().submit(pullFeaturesRunnable).get();
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        long endPersistFeaturesTime = System.currentTimeMillis();

        Long chromosomeDuration = (endPersistChromosomeTime - startPersistChromosomeTime) / 1000;
        logger.info("duration to persist Chromosomes: {} seconds", chromosomeDuration);
        Long genesAndGenomeReferencesDuration = (endPersistGenesAndGenomeReferencesTime - startPersistGenesAndGenomeReferencesTime) / 1000;
        logger.info("duration to persist Genes & GenomeReferences: {} seconds", genesAndGenomeReferencesDuration);
        Long referencesSequencesDuration = (endPersistReferenceSequencesTime - startPersistReferenceSequencesTime) / 1000;
        logger.info("duration to persist ReferenceSequences: {} seconds", referencesSequencesDuration);
        Long alignmentsDuration = (endPersistAlignmentsTime - startPersistAlignmentsTime) / 1000;
        logger.info("duration to persist Alignments: {} seconds", alignmentsDuration);
        Long addAlignmentUTRsDuration = (endAddAlignmentUTRsTime - startAddAlignmentUTRsTime) / 1000;
        logger.info("duration to persist Alignment UTRs: {} seconds", addAlignmentUTRsDuration);
        Long featuresDuration = (endPersistFeaturesTime - startPersistFeaturesTime) / 1000;
        logger.info("duration to persist Features: {} seconds", featuresDuration);

        Long totalDuration = (chromosomeDuration + genesAndGenomeReferencesDuration + referencesSequencesDuration + alignmentsDuration
                + featuresDuration + addAlignmentUTRsDuration) / 60;
        logger.info("Total time to pull from NCBI: {} minutes", totalDuration);

    }

}

package org.renci.hearsay.commands.ncbi;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.renci.gene2accession.G2AFilter;
import org.renci.gene2accession.G2AParser;
import org.renci.gene2accession.filter.G2AAndFilter;
import org.renci.gene2accession.filter.G2AAssemblyFilter;
import org.renci.gene2accession.filter.G2AGenomicNucleotideAccessionVersionPrefixFilter;
import org.renci.gene2accession.filter.G2AProteinAccessionVersionPrefixFilter;
import org.renci.gene2accession.filter.G2ARNANucleotideAccessionVersionPrefixFilter;
import org.renci.gene2accession.filter.G2ATaxonIdFilter;
import org.renci.gene2accession.model.Record;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullReferenceSequencesRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PullReferenceSequencesRunnable.class);

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
            List<G2AFilter> filters = Arrays.asList(new G2AFilter[] { new G2ATaxonIdFilter(9606),
                    new G2AAssemblyFilter("Reference.*(Primary Assembly|ALT_REF_LOCI.*)"),
                    new G2AProteinAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NP_" })),
                    // new G2AGenomicNucleotideAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NC_" })),
                    new G2ARNANucleotideAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })) });
            G2AAndFilter andFilter = new G2AAndFilter(filters);
            List<Record> recordList = gene2AccessionParser.parse(andFilter, genes2RefSeqFile);
            ExecutorService es = Executors.newFixedThreadPool(8);
            for (Record record : recordList) {
                es.submit(new PersistReferenceSequencesRunnable(record, hearsayDAOBean));
            }
            es.shutdown();
            es.awaitTermination(1L, TimeUnit.HOURS);

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

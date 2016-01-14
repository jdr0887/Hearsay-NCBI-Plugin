package org.renci.hearsay.commands.ncbi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.collections4.CollectionUtils;
import org.junit.Test;
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

public class PullAlignmentsTest {

    @Test
    public void testReferenceSequenceRetrieval() {

        // this will take a while
        GBFFManager gbffMgr = GBFFManager.getInstance(1, true);

        List<GBFFFilter> filters = Arrays
                .asList(new GBFFFilter[] { new GBFFSequenceAccessionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })),
                        new GBFFSourceOrganismNameFilter("Homo sapiens"), new GBFFFeatureSourceOrganismNameFilter("Homo sapiens"),
                        new GBFFFeatureTypeNameFilter("CDS"), new GBFFFeatureTypeNameFilter("source") });

        GBFFAndFilter gbffFilter = new GBFFAndFilter(filters);

        List<File> fileList = FTPUtil.ncbiDownloadBySuffix("/refseq/H_sapiens/mRNA_Prot", "rna.gbff.gz");

        fileList.forEach(a -> System.out.println(a.getAbsolutePath()));

        BufferedWriter bw = null;

        try {

            bw = new BufferedWriter(new FileWriter(new File("/tmp", "refseq-protein.txt")));

            for (File f : fileList) {

                System.out.printf("parsing GenBankFlatFile: %s%n", f.getAbsolutePath());

                List<Sequence> sequenceList = gbffMgr.deserialize(gbffFilter, f);

                if (CollectionUtils.isEmpty(sequenceList)) {
                    System.out.println("no sequences found");
                    continue;
                }

                System.out.printf("sequenceList.size(): %s%n", sequenceList.size());

                for (Sequence sequence : sequenceList) {

                    System.out.println(sequence.toString());

                    if (CollectionUtils.isEmpty(sequence.getFeatures())) {
                        System.out.println("sequence.getFeatures() is empty");
                        continue;
                    }

                    // protein accession
                    String proteinAccession = null;
                    Feature firstCDSFeature = null;
                    for (Feature feature : sequence.getFeatures()) {
                        if (!"CDS".equals(feature.getType())) {
                            continue;
                        }
                        firstCDSFeature = feature;
                        break;
                    }
                    proteinAccession = firstCDSFeature.getQualifiers().get("protein_id").replace("\"", "");

                    int exonCount = 0;
                    for (Feature feature : sequence.getFeatures()) {
                        if ("exon".equals(feature.getType())) {
                            exonCount++;
                        }
                    }

                    if (exonCount == 0) {
                        System.out.printf("no exons found: %s%n", sequence.toString());
                        continue;
                    }

                    System.out.printf("number of exons found: %s%n", exonCount);

                    String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ")
                            ? sequence.getVersion().substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();

                    bw.write(String.format("%s|%s%n", proteinAccession, refSeqVersionedAccession));
                    bw.flush();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                bw.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}

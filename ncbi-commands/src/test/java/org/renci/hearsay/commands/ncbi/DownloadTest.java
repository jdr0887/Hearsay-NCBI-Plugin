package org.renci.hearsay.commands.ncbi;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;
import org.renci.gbff.GBFFFilter;
import org.renci.gbff.GBFFManager;
import org.renci.gbff.filter.GBFFAndFilter;
import org.renci.gbff.filter.GBFFFeatureSourceOrganismNameFilter;
import org.renci.gbff.filter.GBFFFeatureTypeNameFilter;
import org.renci.gbff.filter.GBFFSequenceAccessionPrefixFilter;
import org.renci.gbff.filter.GBFFSourceOrganismNameFilter;
import org.renci.gbff.model.Sequence;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;

public class DownloadTest {

    @Test
    public void humanDownload() {

        GBFFManager gbffMgr = GBFFManager.getInstance(1, true);

        List<GBFFFilter> filters = Arrays
                .asList(new GBFFFilter[] { new GBFFSequenceAccessionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })),
                        new GBFFSourceOrganismNameFilter("Homo sapiens"), new GBFFFeatureSourceOrganismNameFilter("Homo sapiens"),
                        new GBFFFeatureTypeNameFilter("CDS"), new GBFFFeatureTypeNameFilter("source") });

        GBFFAndFilter gbffFilter = new GBFFAndFilter(filters);

        List<File> fileList = FTPUtil.ncbiDownloadBySuffix("/refseq/H_sapiens/mRNA_Prot", "rna.gbff.gz");

        List<Sequence> sequenceList = new ArrayList<Sequence>();
        for (File f : fileList) {
            List<Sequence> tmpList = gbffMgr.deserialize(gbffFilter, f);
            System.out.printf("sequences: %d, file: %s%n", tmpList.size(), f.getAbsolutePath());
            sequenceList.addAll(tmpList);
        }
        System.out.println(sequenceList.size());
    }

    @Test
    public void vertebrateMammalianDownload() {

        GBFFManager gbffMgr = GBFFManager.getInstance(8, true);

        List<GBFFFilter> filters = Arrays
                .asList(new GBFFFilter[] { new GBFFSequenceAccessionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })),
                        new GBFFSourceOrganismNameFilter("Homo sapiens"), new GBFFFeatureSourceOrganismNameFilter("Homo sapiens"),
                        new GBFFFeatureTypeNameFilter("CDS"), new GBFFFeatureTypeNameFilter("source") });

        GBFFAndFilter gbffFilter = new GBFFAndFilter(filters);

        List<File> fileList = FTPUtil.ncbiDownloadBySuffix("/refseq/release/vertebrate_mammalian", "rna.gbff.gz");
        List<Sequence> sequenceList = new ArrayList<Sequence>();
        for (File f : fileList) {
            List<Sequence> tmpList = gbffMgr.deserialize(gbffFilter, f);
            System.out.printf("sequences: %d, file: %s%n", tmpList.size(), f.getAbsolutePath());
            sequenceList.addAll(tmpList);
        }
        System.out.println(sequenceList.size());
    }

}

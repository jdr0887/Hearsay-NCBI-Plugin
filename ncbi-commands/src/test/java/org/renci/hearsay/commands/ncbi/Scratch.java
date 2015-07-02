package org.renci.hearsay.commands.ncbi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;

public class Scratch {

    @Test
    public void test() {
        // String asdf = "join(1..3,4..5)";
        String asdf = "order(366..368,375..377,381..383,390..392,399..401)";

        Pattern p = Pattern.compile("^(join|order)\\((.+)\\)$");
        Matcher m = p.matcher(asdf);
        System.out.println(m.find());
        System.out.println(m.group(2));
        // System.out.println(m.group(2));
    }

    @Test
    public void scratch() throws IOException {
        File f = new File("/home/jdr0887/Downloads", "gene2refseq.filtered");
        List<String> lines = FileUtils.readLines(f);
        for (String line : lines) {
            String[] asdf = line.split("\t");
            System.out.println(asdf[0]);
            System.out.println(asdf[9]);
            System.out.println(asdf[10]);
            System.out.println(asdf[12]);
            System.out.println(asdf[15]);
        }
    }

    @Test
    public void testDownload() {
        File genes2RefSeqFile = FTPUtil.ncbiDownload("/gene/DATA", "gene2refseq.gz");
        Map<String, Pair<Integer, Integer>> genes2RefSeqMap = new HashMap<String, Pair<Integer, Integer>>();
        // #Format: tax_id GeneID status RNA_nucleotide_accession.version RNA_nucleotide_gi protein_accession.version
        // protein_gi genomic_nucleotide_accession.version genomic_nucleotide_gi start_position_on_the_genomic_accession
        // end_position_on_the_genomic_accession orientation assembly mature_peptide_accession.version mature_peptide_gi
        // Symbol

        try (FileInputStream fis = new FileInputStream(genes2RefSeqFile);
                GZIPInputStream gis = new GZIPInputStream(fis);
                InputStreamReader isr = new InputStreamReader(gis);
                BufferedReader br = new BufferedReader(isr)) {
            String line;
            while ((line = br.readLine()) != null) {

                if (line.startsWith("#")) {
                    continue;
                }

                Scanner scanner = new Scanner(line).useDelimiter("\t");
                // String[] split = line.split("\t");

                // tax_id
                String taxId = scanner.next();
                // protein_gi
                scanner.next();
                // genomic_nucleotide_accession.version
                scanner.next();
                // GeneID
                scanner.next();
                // status
                scanner.next();
                // RNA_nucleotide_accession.version
                String rnaNucleotideAccession = scanner.next();
                // RNA_nucleotide_gi
                scanner.next();
                // protein_accession.version
                scanner.next();
                // genomic_nucleotide_gi
                scanner.next();
                // start_position_on_the_genomic_accession
                String genomicStart = scanner.next();
                // end_position_on_the_genomic_accession
                String genomicEnd = scanner.next();
                // orientation
                scanner.next();
                // assembly mature_peptide_accession.version mature_peptide_gi Symbol
                String assembly = scanner.next();
                // mature_peptide_accession.version
                scanner.next();
                // mature_peptide_gi
                scanner.next();
                // Symbol
                scanner.next();

                // only want human
                if (!"9606".equals(taxId)) {
                    continue;
                }

                if (!assembly.contains("GRCh38.p2")) {
                    continue;
                }

                if (!assembly.contains("Primary Assembly")) {
                    continue;
                }

                if (!genes2RefSeqMap.containsKey(rnaNucleotideAccession)) {
                    genes2RefSeqMap.put(rnaNucleotideAccession,
                            new Pair<Integer, Integer>(Integer.valueOf(genomicStart), Integer.valueOf(genomicEnd)));
                }
                scanner.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}

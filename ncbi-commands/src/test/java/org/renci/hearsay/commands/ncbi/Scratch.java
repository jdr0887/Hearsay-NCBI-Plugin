package org.renci.hearsay.commands.ncbi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;
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
import org.renci.hearsay.dao.model.DirectionType;

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
        Pattern p = Pattern.compile("(\\d+)([-|+])(\\d+)");
        Matcher m = p.matcher("127+12");
        if (m.matches()) {
            DirectionType directionType = null;
            for (DirectionType dt : DirectionType.values()) {
                if (dt.getValue().equals(m.group(2))) {
                    directionType = dt;
                    break;
                }
            }
            System.out.println(directionType.getValue());
            System.out.println(m.group(1));
            System.out.println(m.group(3));
        }
    }

    @Test
    public void testGene2RefSeqFilters() {

        File genes2RefSeqFile = FTPUtil.ncbiDownload("/gene/DATA", "gene2refseq.gz");
        G2AParser gene2AccessionParser = G2AParser.getInstance(8);
        List<G2AFilter> filters = Arrays.asList(new G2AFilter[] { new G2ATaxonIdFilter(9606),
                // new G2AAssemblyFilter("Reference.*(Primary Assembly|ALT_REF_LOCI.*)"),
                new G2AAssemblyFilter("Reference.*Primary Assembly"),
                new G2AProteinAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NP_" })),
                new G2AGenomicNucleotideAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NC_" })),
                new G2ARNANucleotideAccessionVersionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })) });
        G2AAndFilter andFilter = new G2AAndFilter(filters);
        List<Record> recordList = gene2AccessionParser.parse(andFilter, genes2RefSeqFile);
        System.out.println(recordList.size());

    }

    @Test
    public void testRefseqAssemblySummaryFile() throws IOException {

        try (BufferedReader br = new BufferedReader(new FileReader(new File("/tmp", "assembly_summary_historical.txt")))) {
            // # assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid
            // organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep
            // seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                }
                try (Scanner scanner = new Scanner(line).useDelimiter("\t")) {
                    String assemblyAccession = scanner.next();
                    String bioProject = scanner.next();
                    String bioSample = scanner.next();
                    String wgsMaster = scanner.next();
                    String refseqCategory = scanner.next();
                    String taxId = scanner.next();
                    String speciesTaxId = scanner.next();
                    String organismName = scanner.next();
                    if (!"homo sapiens".equalsIgnoreCase(organismName)) {
                        continue;
                    }
                    String infraspecificName = scanner.next();
                    String isolate = scanner.next();
                    String versionStatus = scanner.next();
                    String assemblyLevel = scanner.next();
                    String releaseType = scanner.next();
                    String genomeRep = scanner.next();
                    String seqReleaseDate = scanner.next();
                    String asmName = scanner.next();
                    if (!asmName.startsWith("GR")) {
                        continue;
                    }
                    String submitter = scanner.next();
                    String gbrsPairedASM = scanner.next();
                    String pairedASMComp = scanner.next();
                    String ftpPath = scanner.next();
                    System.out.println(assemblyAccession);
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

package org.renci.hearsay.commands.ncbi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.junit.Test;

public class PullGenesTest {

    @Test
    public void test() {
        File genesFile = new File("/home/jdr0887/Downloads", "Homo_sapiens.gene_info.gz");

        Set<String> dataSet = new HashSet<String>();

        // parse
        try (FileInputStream fis = new FileInputStream(genesFile);
                GZIPInputStream gis = new GZIPInputStream(fis);
                InputStreamReader isr = new InputStreamReader(gis);
                BufferedReader br = new BufferedReader(isr)) {

            // #Format: tax_id GeneID Symbol LocusTag Synonyms dbXrefs chromosome map_location description type_of_gene
            // Symbol_from_nomenclature_authority Full_name_from_nomenclature_authority Nomenclature_status
            // Other_designations Modification_date (tab is used as a separator, pound sign - start of a comment)
            String line;
            while ((line = br.readLine()) != null) {

                if (line.startsWith("#")) {
                    continue;
                }

                try (Scanner scanner = new Scanner(line).useDelimiter("\t")) {

                    String taxId = scanner.next();
                    String geneId = scanner.next();
                    String symbol = scanner.next();
                    String locusTag = scanner.next();
                    String synonyms = scanner.next();
                    String dbXrefs = scanner.next();
                    String chromosome = scanner.next();
                    String mapLocation = scanner.next();
                    String description = scanner.next();
                    // String typeOfGene = st.nextToken();
                    // String symbolFromNomenclatureAuthority = st.nextToken();
                    // String nameFromNomenclatureAuthority = st.nextToken();
                    // String nomenclatureStatus = st.nextToken();
                    // String otherDesignations = st.nextToken();
                    // String modificationDate = st.nextToken();

                    if (chromosome.equals("-") || chromosome.equalsIgnoreCase("Un") || !chromosome.contains("|")) {
                        continue;
                    }

                    dataSet.add(String.format("%s\t%s", symbol, chromosome));

                }

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        for (String data : dataSet) {
            System.out.println(data);
        }

    }

}

package org.renci.hearsay.commands.ncbi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

import org.apache.commons.collections4.CollectionUtils;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.model.Chromosome;
import org.renci.hearsay.dao.model.Gene;
import org.renci.hearsay.dao.model.GeneSymbol;
import org.renci.hearsay.dao.model.Identifier;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullGenesRunnable implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(PullGenesRunnable.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullGenesRunnable() {
        super();
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        File genesFile = FTPUtil.ncbiDownload("/gene/DATA/GENE_INFO/Mammalia", "Homo_sapiens.gene_info.gz");

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

                    if (chromosome.equals("-") || chromosome.equalsIgnoreCase("Un")) {
                        continue;
                    }

                    Gene gene = new Gene();
                    gene.setDescription(description);
                    gene.setSymbol(symbol);
                    List<Gene> potentiallyFoundGeneList = hearsayDAOBean.getGeneDAO().findByExample(gene);
                    if (CollectionUtils.isNotEmpty(potentiallyFoundGeneList)) {
                        logger.warn("Gene is already persisted: {}", symbol);
                        continue;
                    }
                    gene.setId(hearsayDAOBean.getGeneDAO().save(gene));
                    logger.info(gene.toString());

                    if (chromosome.indexOf("|") != -1) {
                        String[] split = chromosome.split("|");
                        for (String chr : split) {
                            List<Chromosome> potentialChromosomeList = hearsayDAOBean.getChromosomeDAO()
                                    .findByName(chr);
                            if (CollectionUtils.isNotEmpty(potentialChromosomeList)) {
                                gene.getChromosomes().addAll(potentialChromosomeList);
                            }
                        }
                    } else {
                        List<Chromosome> potentialChromosomeList = hearsayDAOBean.getChromosomeDAO().findByName(
                                chromosome);
                        if (CollectionUtils.isNotEmpty(potentialChromosomeList)) {
                            gene.getChromosomes().addAll(potentialChromosomeList);
                        }
                    }

                    hearsayDAOBean.getGeneDAO().save(gene);

                    if (!synonyms.trim().equals("-")) {
                        StringTokenizer geneSymbolStringTokenizer = new StringTokenizer(synonyms, "|");

                        while (geneSymbolStringTokenizer.hasMoreTokens()) {
                            String geneSymbol = geneSymbolStringTokenizer.nextToken();
                            GeneSymbol gs = new GeneSymbol();
                            gs.setSymbol(geneSymbol);
                            gs.setGene(gene);
                            gs.setId(hearsayDAOBean.getGeneSymbolDAO().save(gs));
                            logger.debug(geneSymbol.toString());
                            gene.getAliases().add(gs);
                        }
                    }

                    Identifier identifier = new Identifier("www.ncbi.nlm.nih.gov/gene", geneId);
                    identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                    logger.debug(identifier.toString());
                    gene.getIdentifiers().add(identifier);
                    hearsayDAOBean.getGeneDAO().save(gene);

                }

            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        // genesFile.delete();
    }

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}

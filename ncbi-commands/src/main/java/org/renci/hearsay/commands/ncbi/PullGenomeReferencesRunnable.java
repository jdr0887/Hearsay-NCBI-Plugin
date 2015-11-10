package org.renci.hearsay.commands.ncbi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.List;
import java.util.Scanner;

import org.apache.commons.collections4.CollectionUtils;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.HearsayDAOException;
import org.renci.hearsay.dao.model.GenomeReference;
import org.renci.hearsay.dao.model.Identifier;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullGenomeReferencesRunnable implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(PullGenomeReferencesRunnable.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullGenomeReferencesRunnable() {
        super();
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        File refseqAssemblySummaryFile = FTPUtil.ncbiDownload("/genomes/refseq", "assembly_summary_refseq.txt");

        try (BufferedReader br = new BufferedReader(new FileReader(refseqAssemblySummaryFile))) {
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

                    GenomeReference exampleGenomeReference = new GenomeReference();
                    exampleGenomeReference.setName(asmName);

                    List<GenomeReference> potentiallyFoundGenomeReferenceList = hearsayDAOBean.getGenomeReferenceDAO()
                            .findByExample(exampleGenomeReference);
                    logger.info(exampleGenomeReference.toString());
                    if (CollectionUtils.isNotEmpty(potentiallyFoundGenomeReferenceList)) {
                        logger.info("GenomeReference is already persisted");
                        continue;
                    }

                    exampleGenomeReference.setId(hearsayDAOBean.getGenomeReferenceDAO().save(exampleGenomeReference));

                    String id = getGenomeReferenceAssemblyId(assemblyAccession);

                    Identifier identifier = new Identifier("www.ncbi.nlm.nih.gov/assembly", id);
                    identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                    logger.info(identifier.toString());

                    exampleGenomeReference.getIdentifiers().add(identifier);

                    hearsayDAOBean.getGenomeReferenceDAO().save(exampleGenomeReference);

                }
            }
        } catch (HearsayDAOException | IOException e) {
            e.printStackTrace();
        }
        // refseqAssemblySummaryFile.delete();
    }

    private String getGenomeReferenceAssemblyId(String assemblyAccession) {
        String ret = null;
        try {
            String url = String.format("http://www.ncbi.nlm.nih.gov/assembly/%s?report=xml&format=text",
                    assemblyAccession);
            URL obj = new URL(url);
            HttpURLConnection con = (HttpURLConnection) obj.openConnection();
            con.setRequestMethod("GET");
            con.setRequestProperty("User-Agent", "Mozilla/5.0");

            int responseCode = con.getResponseCode();
            logger.info("\nSending 'GET' request to URL : " + url);
            logger.info("Response Code : " + responseCode);

            BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));
            String inputLine;
            while ((inputLine = in.readLine()) != null) {
                if (inputLine.contains("DocumentSummary")) {
                    ret = inputLine.substring(inputLine.indexOf("\"") + 1, inputLine.lastIndexOf("\""));
                    break;
                }
            }
            in.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return ret;
    }

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}

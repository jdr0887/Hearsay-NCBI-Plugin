package org.renci.hearsay.commands.clinvar;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.List;

import junit.framework.Assert;

import org.apache.commons.io.FileUtils;
import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import org.apache.commons.net.ftp.FTPReply;
import org.junit.Test;

public class Scratch {

    @Test
    public void getVersion() {

        String ncbiHost = "ftp.ncbi.nlm.nih.gov";
        FTPClient ftpClient = new FTPClient();
        try {

            ftpClient.connect(ncbiHost);
            ftpClient.login("anonymous", "anonymous");
            ftpClient.setFileType(FTP.ASCII_FILE_TYPE);
            ftpClient.enterLocalPassiveMode();
            int reply = ftpClient.getReplyCode();
            if (!FTPReply.isPositiveCompletion(reply)) {
                ftpClient.disconnect();
                System.err.println("FTP server refused connection.");
            }

            File tmpFile = new File(System.getProperty("java.io.tmpdir", "/tmp"), "version.txt");
            OutputStream fos = new BufferedOutputStream(new FileOutputStream(tmpFile));
            ftpClient.retrieveFile("/genomes/H_sapiens/README_CURRENT_RELEASE", fos);
            fos.flush();
            fos.close();

            String version = "";
            List<String> lines = FileUtils.readLines(tmpFile);
            for (String line : lines) {
                if (line.startsWith("ASSEMBLY NAME:")) {
                    version = line.substring(line.indexOf(":") + 1, line.length()).trim();
                    break;
                }
            }

            System.out.println("version: " + version);

            List<FTPFile> files = Arrays.asList(ftpClient.listFiles("/genomes/H_sapiens/Assembled_chromosomes/seq/"));

            for (FTPFile file : files) {
                if (file.getName().startsWith(String.format("hs_ref_%s", version)) && file.getName().endsWith(".fa.gz")) {
                    File f = new File(System.getProperty("java.io.tmpdir", "/tmp"), file.getName());
                    fos = new BufferedOutputStream(new FileOutputStream(f));
                    ftpClient.retrieveFile(
                            String.format("/genomes/H_sapiens/Assembled_chromosomes/seq/%s", file.getName()), fos);
                    fos.flush();
                    fos.close();
                }
            }

        } catch (IOException e) {

            if (ftpClient.isConnected()) {
                try {
                    ftpClient.disconnect();
                } catch (IOException f) {
                    // do nothing
                }
            }
            System.err.println("Could not connect to server.");
            e.printStackTrace();
        }

    }

}

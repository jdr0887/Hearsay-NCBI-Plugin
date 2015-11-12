package org.renci.hearsay.commands.ncbi;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;

import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPReply;
import org.renci.lrg.Lrg;

public class PullLocusReferenceGenomicAction {

    public PullLocusReferenceGenomicAction() {
        super();
    }

    public void pull() {

        File tmpFile = new File(System.getProperty("java.io.tmpdir", "/tmp"), "LRG_public_xml_files.zip");
        FTPClient ftpClient = new FTPClient();
        // download
        try {
            ftpClient.connect("ftp.ebi.ac.uk");

            ftpClient.login("anonymous", "anonymous");
            ftpClient.setFileType(FTP.BINARY_FILE_TYPE);
            ftpClient.enterLocalPassiveMode();

            int reply = ftpClient.getReplyCode();
            if (!FTPReply.isPositiveCompletion(reply)) {
                ftpClient.disconnect();
                System.err.println("FTP server refused connection.");
                return;
            }

            try (OutputStream fos = new BufferedOutputStream(new FileOutputStream(tmpFile))) {
                ftpClient.retrieveFile("/pub/databases/lrgex/LRG_public_xml_files.zip", fos);
                fos.flush();
            }

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (ftpClient.isConnected()) {
                    ftpClient.disconnect();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        // parse
        try (ZipFile zipFile = new ZipFile(tmpFile)) {

            final Enumeration<? extends ZipEntry> entries = zipFile.entries();
            while (entries.hasMoreElements()) {
                final ZipEntry entry = entries.nextElement();
                System.out.println(entry.getName());
                // use entry input stream:
                try (InputStream is = zipFile.getInputStream(entry)) {
                    JAXBContext context = JAXBContext.newInstance(Lrg.class);
                    Unmarshaller unmarshaller = context.createUnmarshaller();
                    Lrg lrg = (Lrg) unmarshaller.unmarshal(is);
                    System.out.println(lrg.getFixedAnnotation().getId());
                } catch (JAXBException e) {
                    e.printStackTrace();
                }
            }

        } catch (IOException e1) {
            e1.printStackTrace();
        }
    }

}

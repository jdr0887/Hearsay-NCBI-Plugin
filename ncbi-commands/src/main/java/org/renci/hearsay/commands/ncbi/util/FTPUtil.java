package org.renci.hearsay.commands.ncbi.util;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPReply;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class FTPUtil {

    private static final Logger logger = LoggerFactory.getLogger(FTPUtil.class);

    public static File download(String host, String path, String name) {
        logger.info("downloading: {}", String.format("%s:%s/%s", host, path, name));
        File ret = new File("/tmp", name);
        if (ret.exists()) {
            return ret;
        }
        ret = new File(System.getProperty("java.io.tmpdir", "/tmp"), name);
        FTPClient ftpClient = new FTPClient();
        try {
            ftpClient.connect(host);
            ftpClient.login("anonymous", "anonymous");
            ftpClient.setFileType(FTP.BINARY_FILE_TYPE);
            ftpClient.enterLocalPassiveMode();
            int reply = ftpClient.getReplyCode();
            if (!FTPReply.isPositiveCompletion(reply)) {
                ftpClient.disconnect();
                logger.error("FTP server refused connection.");
                return null;
            }
            try (OutputStream fos = new BufferedOutputStream(new FileOutputStream(ret))) {
                ftpClient.retrieveFile(String.format("%s/%s", path, name), fos);
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
        return ret;
    }

    public static File ncbiDownload(String path, String name) {
        String host = "ftp.ncbi.nlm.nih.gov";
        return download(host, path, name);
    }

    public static File ucscDownload(String path, String name) {
        String host = "hgdownload.cse.ucsc.edu";
        return download(host, path, name);
    }

}

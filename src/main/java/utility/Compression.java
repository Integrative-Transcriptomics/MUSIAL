package utility;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.zip.GZIPOutputStream;

/**
 * Implements static methods to compress {@link String} content using the {@link java.util.zip} library.
 *
 * @author Simon Hackl
 */
public class Compression {

    /**
     * Compress the byte content of a {@link String} using `gzip` compression and returns the bytes of the compressed data.
     *
     * @param content {@link String} to compress.
     * @return {@link Byte} Array; `gzip` compressed byte content of the input.
     * @throws IOException See {@link GZIPOutputStream#write} and {@link GZIPOutputStream#close()}.
     */
    public static byte[] gzip(String content) throws IOException {
        ByteArrayOutputStream outputStream = new ByteArrayOutputStream(512);
        GZIPOutputStream gzipOutputStream = new GZIPOutputStream(outputStream);
        gzipOutputStream.write(content.getBytes());
        gzipOutputStream.close();
        return outputStream.toByteArray();
    }

}

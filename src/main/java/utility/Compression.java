package utility;

import com.aayushatharva.brotli4j.Brotli4jLoader;
import com.aayushatharva.brotli4j.decoder.Decoder;
import com.aayushatharva.brotli4j.encoder.Encoder;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Base64;

/**
 * Comprises static methods to compress {@link String} content using the `brotli4j` library.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.2
 */
public class Compression {

    /**
     * Compress and encode the byte content of a {@link String} using `brotli` compression and returns the Base64 encoded {@link String}
     * representation of the compressed data.
     *
     * @param content {@link String} to compress.
     * @return {@link String}; Base64 representation of the `brotli` compressed byte content of the input.
     * @throws IOException .
     */
    public static String brotliEncodeString(String content) throws IOException {
        Brotli4jLoader.ensureAvailability();
        byte[] byteContent = content.getBytes();
        byte[] byteContentCompressed = Encoder.compress(byteContent);
        return Base64.getEncoder().encodeToString(byteContentCompressed);
    }

    /**
     * Compress and encode the byte content of a {@link String} using `brotli` compression and returns the byte content of the compressed data.
     *
     * @param content {@link String} to compress.
     * @return {@link Byte} Array; `brotli` compressed byte content of the input.
     * @throws IOException .
     */
    public static byte[] brotliEncodeStringToBytes(String content) throws IOException {
        Brotli4jLoader.ensureAvailability();
        byte[] byteContent = content.getBytes();
        return Encoder.compress(byteContent);
    }

    /**
     * Decompress and decode the `brotli` compressed Base64-representation byte content of a {@link String}.
     *
     * @param content {@link String} to decode.
     * @return {@link String}; Decoded content of the `brotli` compressed input.
     * @throws IOException .
     */
    public static String brotliDecodeString(String content) throws IOException {
        return brotliDecodeBytes(Base64.getDecoder().decode(content));
    }

    /**
     * Decompress and decode `brotli` compressed byte content.
     *
     * @param content {@link Byte} array to decode.
     * @return {@link String}; Decoded content of the `brotli` compressed input.
     * @throws IOException .
     */
    public static String brotliDecodeBytes(byte[] content) throws IOException {
        Brotli4jLoader.ensureAvailability();
        byte[] byteContentDecompressed = Decoder.decompress(content).getDecompressedData();
        return new String(byteContentDecompressed, StandardCharsets.UTF_8);
    }

}

package op;

import com.google.common.base.Splitter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import main.Musial;
import model.Contig;
import model.Feature;
import model.Storage;
import util.Bio;
import util.Constants;
import util.Logging;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;

/**
 * The {@code StorageIO} class provides utility methods for serializing and deserializing genomic data.
 * <p>
 * This class includes methods to convert {@link Storage} objects into various file formats such as JSON, GFF3, FASTA, and VCF. It handles
 * the generation of file content based on the data stored in the {@link Storage} object, ensuring compliance with the respective file
 * format specifications. Additionally, it provides helper methods for processing features and contigs.
 */
public class StorageIO {

    /**
     * Private constructor to prevent instantiation of this utility class.
     */
    private StorageIO() {
    }

    /**
     * Gson instance for JSON serialization and deserialization with pretty printing.
     */
    private static final Gson gson = new GsonBuilder()
            .setPrettyPrinting()
            .create();

    /**
     * Serializes the given {@link Storage} object to a JSON file at the specified path.
     * <p>
     * This method converts the {@link Storage} object into a JSON string using the Gson library. If the specified path does not end with
     * ".json" or ".json.gz", the default output extension defined in {@link Musial} is appended to the path. The JSON data is then written
     * to the file.
     * <p>
     * If the path ends with ".gz", the JSON data is compressed using GZIP before being written. Otherwise, it is written as plain text.
     * </p>
     *
     * @param storage The {@link Storage} object to be serialized.
     * @param path    The {@link Path} where the JSON file will be written.
     * @throws IOException If an I/O error occurs during file writing.
     */
    public static void toJSON(Storage storage, Path path) throws IOException {
        // Ensure the file path ends with a valid extension.
        if (!(path.toString().endsWith(".json") || path.toString().endsWith(".json.gz"))) {
            path = Path.of(path + Musial.OUTPUT_EXTENSION);
        }

        // Serialize the Storage object to a JSON string.
        String jsonData = gson.toJson(storage, Storage.class);

        // Write the JSON data to the specified file, compressing it if necessary.
        try (Writer writer = path.toString().endsWith(".gz")
                ? new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(path.toFile())))
                : new FileWriter(path.toFile(), StandardCharsets.UTF_8)) {
            writer.write(jsonData);
        } catch (IOException e) {
            // Throw an exception with a detailed error message if writing fails.
            throw new IOException(String.format("Failed to write MUSIAL storage to file %s; %s.", path,
                    e.getMessage()));
        }
    }

    /**
     * Generates the content of a GFF (General Feature Format) file from the given {@link Storage} object.
     * <p>
     * This method constructs a GFF file content as a {@link String} by iterating over the features in the provided {@link Storage} object.
     * The GFF content includes the version, processor information, and the feature data. Each feature is converted to its GFF string
     * representation using the {@link #featureToGFF3String(Feature)} method.
     * <p>
     * The generated GFF content follows the GFF3 specification and includes the following:
     * <ul>
     *   <li>##gff-version: Specifies the GFF version.</li>
     *   <li>##processor: Includes the software id and version used to generate the file.</li>
     *   <li>Feature data: Each feature is represented in GFF format.</li>
     * </ul>
     *
     * @param storage The {@link Storage} object containing the features to include in the GFF file.
     * @return A {@link String} representing the GFF file content.
     */
    public static String toGFF3(Storage storage) {
        StringBuilder content = new StringBuilder();
        content.append("##gff-version 3.1.26").append(Constants.LINE_SEPARATOR);
        content.append("##processor %s %s".formatted(Musial.name, Musial.version)).append(Constants.LINE_SEPARATOR);
        for (Feature feature : storage.getFeatures()) {
            content.append(featureToGFF3String(feature));
        }
        return content.toString();
    }

    /**
     * Generates the content of a FASTA file from the given {@link Storage} object.
     * <p>
     * This method constructs a FASTA file content as a {@link String} by iterating over the contigs in the provided {@link Storage} object.
     * Each contig's ID is used as the header (prefixed with '>'), and its sequence is split into lines of 80 characters for proper FASTA
     * formatting. The method ensures that all contigs in the storage have sequence data before proceeding.
     * </p>
     *
     * @param storage The {@link Storage} object containing the contigs and their sequences.
     * @return A {@link String} representing the content of the reference FASTA file.
     * @throws IOException              If an I/O error occurs during the generation of the FASTA content.
     * @throws IllegalArgumentException If no reference sequence information is stored in the {@link Storage} object.
     */
    public static String toFASTA(Storage storage) throws IOException {
        // Validate that all contigs in the storage have sequence data.
        if (!storage.getContigs().stream().allMatch(Contig::hasSequence)) {
            throw new IllegalArgumentException("No reference sequence information is stored in the specified storage.");
        }

        // Initialize a StringBuilder to construct the FASTA content.
        StringBuilder content = new StringBuilder();

        // Iterate over each contig in the storage.
        for (Contig contig : storage.getContigs()) {
            // Append the contig ID as the FASTA header, prefixed with '>'.
            content.append(">").append(contig._id).append(Constants.LINE_SEPARATOR);

            // Split the contig sequence into lines of 80 characters and append each line.
            Splitter.fixedLength(80).split(contig.getSequence()).forEach(line ->
                    content.append(line).append(Constants.LINE_SEPARATOR)
            );
        }

        // Return the constructed FASTA content as a string.
        return content.toString();
    }

    /**
     * Generates the content of a VCF (Variant Call Format) file from the given {@link Storage} object.
     * <p>
     * This method constructs a VCF file content as a {@link String} by iterating over the contigs in the provided {@link Storage} object.
     * The VCF content includes the file format, source, and a header line, followed by the variant data. Each variant is represented by its
     * chromosome, position, reference base, and alternate base.
     * <p>
     * The generated VCF content follows the VCFv4.3 specification and includes the following fields:
     * <ul>
     *   <li>CHROM: Chromosome identifier.</li>
     *   <li>POS: Position of the variant on the chromosome.</li>
     *   <li>ID: Variant identifier (set to ".").</li>
     *   <li>REF: Reference base(s) (gaps are stripped).</li>
     *   <li>ALT: Alternate base(s) (gaps are stripped).</li>
     *   <li>QUAL: Quality score (set to "100").</li>
     *   <li>FILTER: Filter status (set to ".").</li>
     *   <li>INFO: Additional information (see parameters).</li>
     * </ul>
     * <p>
     * Variants can be filtered based on their novelty and ambiguity:
     * <ul>
     *   <li>If {@code onlyNovel} is {@code true}, only active variants are included.</li>
     *   <li>If {@code excludeAmbiguous} is {@code true}, variants with ambiguous alternate bases are excluded.</li>
     * </ul>
     *
     * @param storage                The {@link Storage} object containing the contigs and variants.
     * @param onlyNovel              If {@code true}, only active variants are included in the VCF content.
     * @param excludeAmbiguous       If {@code true}, variants with ambiguous alternate bases are excluded.
     * @return A {@link String} representing the VCF file content.
     */
    public static String toVCF(Storage storage, boolean onlyNovel, boolean excludeAmbiguous) {
        // Initialize the VCF content with the file format, source, and header lines.
        StringBuilder content = new StringBuilder()
                .append("##fileformat=VCFv4.3").append(Constants.LINE_SEPARATOR)
                .append("##source=MUSIAL").append(Constants.LINE_SEPARATOR)
                .append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").append(Constants.LINE_SEPARATOR);

        // Iterate over each contig in the storage.
        storage.getContigs().forEach(contig -> {
            // Retrieve the list of variants based on the novelty filter.
            var variants = onlyNovel ? contig.getActiveVariants() : contig.getAllVariants();

            // Filter and process each variant.
            variants.stream()
                    .filter(variant -> !(excludeAmbiguous && variant.alternative.equals(Constants.ANY_NUCLEOTIDE))) // Exclude ambiguous
                    // variants if required.
                    .forEach(variant -> content.append(String.join("\t",
                                    contig._id, // Chromosome identifier.
                                    String.valueOf(variant.position), // Variant position.
                                    ".", // Variant ID (set to ".").
                                    Bio.stripGaps(variant.reference), // Reference base(s) with gaps stripped.
                                    Bio.stripGaps(variant.alternative), // Alternate base(s) with gaps stripped.
                                    "100", // Quality score.
                                    ".", // Filter status.
                                    "ALT=%s".formatted(variant.alternative))) // Additional information (empty).
                            .append(Constants.LINE_SEPARATOR)); // Append a new line for each variant.
        });

        // Return the constructed VCF content as a string.
        return content.toString();
    }

    /**
     * Converts a {@link Feature} object into its GFF3 (General Feature Format) string representation.
     * <p>
     * This method generates a GFF3-compliant string for the given {@link Feature}, including its attributes and sub-features. The generated
     * string contains the feature's contig, type, start and end positions, strand, and additional attributes. If the feature has
     * sub-features, they are appended to the output.
     * </p>
     *
     * @param feature The {@link Feature} object to be converted to a GFF3 string.
     * @return A {@link String} representing the GFF3 format of the feature and its sub-features.
     */
    private static String featureToGFF3String(Feature feature) {
        StringBuilder contentBuilder = new StringBuilder();

        // Determine GFF3 conforming ID attribute.
        String id = feature.type.contains("gene")
                ? "gene-%s".formatted(feature._id)
                : "%s-%s".formatted(feature.type, feature._id);

        // Log a warning if the feature type is not a gene.
        if (!feature.type.contains("gene")) {
            Logging.logWarning("Feature %s type is not a gene; This may conflict with the GFF3 definition.".formatted(feature.name));
        }

        // Append the main feature information in GFF3 format.
        contentBuilder.append(String.join(Constants.TAB,
                feature.contig, Musial.name, feature.type,
                String.valueOf(feature.start), String.valueOf(feature.end),
                Constants.DOT, String.valueOf(feature.strand), Constants.DOT,
                "ID=%s".formatted(id)));

        // Append feature attributes if present.
        feature.getAttributes().forEach((key, value) ->
                contentBuilder.append(";%s=%s".formatted(key, value)));
        contentBuilder.append(Constants.LINE_SEPARATOR);

        // Append sub-feature information in GFF3 format.
        feature.getSubFeatures().forEach(subFeature ->
                contentBuilder.append(String.join(Constants.TAB,
                                feature.contig, Musial.name, subFeature.type(),
                                String.valueOf(subFeature.start()), String.valueOf(subFeature.end()),
                                Constants.DOT, String.valueOf(feature.strand), Constants.DOT,
                                subFeature.ID(feature._id)))
                        .append(Constants.LINE_SEPARATOR));

        return contentBuilder.toString();
    }

}

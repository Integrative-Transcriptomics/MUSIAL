package datastructure;

import exceptions.MusialException;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Internal representation of a sample, i.e., a set of variant calls from a single biological sample.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.0
 */
@SuppressWarnings("unused")
public final class Sample {

    /**
     * The internal name of the sample.
     */
    public final String name;
    /**
     * A {@link HashMap} yielding any meta-information {@link String} key-value pairs about this {@link Sample}.
     */
    private final HashMap<String, String> annotations = new HashMap<>();
    /**
     * The vcf file that is associated with this sample from which all variant calls are parsed.
     */
    public transient final File vcfFile;
    /**
     * {@link VCFFileReader} instance pointing to the vcf file of this sample.
     */
    public transient VCFFileReader vcfFileReader;

    /**
     * Constructor of {@link Sample}.
     *
     * @param vcfFile {@link File} instances pointing to the .vcf file of this sample.
     * @param name    {@link String} representing the sample name.
     */
    public Sample(File vcfFile, String name) {
        this.vcfFile = vcfFile;
        this.name = name;
    }

    /**
     * Initialize a {@link VCFFileReader} instance pointing to the .vcf file associated with this instance.
     * <p>
     * Validate that a .tbi index of the .vcf file is present and generates one, if not.
     */
    public void imputeVcfFileReader() throws MusialException {
        try {
            // Check for the existence of a `.tbi.gz` index file of the input `.vcf` file.
            if (!new File(this.vcfFile.getAbsolutePath() + ".tbi").exists()) {
                // If none is present, an index is created and written to the same directory as the input `.vcf` file.
                TabixIndex tabixIndex = IndexFactory.createTabixIndex(
                        this.vcfFile,
                        new VCFCodec(),
                        null
                );
                tabixIndex.write(Path.of(this.vcfFile.getAbsolutePath() + ".tbi"));
            }
            // VCFFileReader can now be initialized.
            this.vcfFileReader = new VCFFileReader(
                    this.vcfFile,
                    new File(this.vcfFile.getAbsolutePath() + ".tbi")
            );
        } catch (Exception e) {
            throw new MusialException("Failed to initialize .vcf file reader for sample " + this.name + "; " + e.getMessage());
        }
    }

    /**
     * Closes and removes the {@link VCFFileReader} of this instance, if imputed.
     */
    public void killVcfFileReader() {
        if (this.vcfFileReader != null) {
            this.vcfFileReader.close();
            this.vcfFileReader = null;
        }
    }

    /**
     * Add an annotation to this {@link #annotations}.
     *
     * @param key   Key of the annotation to add.
     * @param value Value of the annotation to add.
     */
    public void addAnnotation(String key, String value) {
        this.annotations.put(key, value);
    }

    /**
     * Remove an annotation from {@link #annotations}, if it exists.
     *
     * @param key Key of the annotation to remove.
     */
    public void removeAnnotation(String key) {
        this.annotations.remove(key);
    }

    /**
     * Query an annotation value from {@link #annotations} or null, if it does not exist.
     *
     * @param key Key of the annotation to query.
     * @return Value of the annotation to query or null.
     */
    public String getAnnotation(String key) {
        return this.annotations.get(key);
    }

    /**
     * Queries all annotation values from {@link #annotations}.
     *
     * @return Set of key/value pairs; One for each annotation stored in {@link #annotations}.
     */
    public Set<Map.Entry<String, String>> getAnnotations() {
        return this.annotations.entrySet();
    }

}

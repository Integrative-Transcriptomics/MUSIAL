package datastructure;

import components.Logging;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;

/**
 * Internal representation of a sample that is subject to analysis.
 * <p>
 * Pairs of instances of this class and {@link FeatureEntry} instances are used internally as so called 'Run
 * entries' to specify a set of analysis tasks, i.e. which sample is analyzed with respect to which specified
 * reference feature.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class SampleEntry {

    /**
     * The `.vcf` format file that is associated with this sample from which all variant information is parsed.
     */
    public transient final File vcfFile;
    /**
     * A {@link VCFFileReader} instance pointing to the `.vcf` file of the respective sample.
     */
    public transient VCFFileReader vcfFileReader;
    /**
     * The internal name of the sample.
     */
    public final String name;
    /**
     * Map of {@link String} key-value pairs; Annotations of this sample. Especially contains one annotation per feature
     * that specifies this samples genotype and proteoform, if the feature is a coding sequence.
     */
    public final HashMap<String, String> annotations = new HashMap<>();
    public final static String REFERENCE_SAMPLE_ID = "WILDTYPE";

    /**
     * Constructor of {@link SampleEntry}.
     *
     * @param vcfReader {@link VCFFileReader} instances initialized with the `.vcf` file of the respective sample.
     * @param name      {@link String} representing the sample name.
     */
    public SampleEntry(File vcfReader, String name) {
        this.vcfFile = vcfReader;
        this.name = name;
    }

    /**
     * Initializes a {@link VCFFileReader} instance pointing to the .vcf file of the specified {@link SampleEntry}.
     * <p>
     * Validates that a Tabix index of the .vcf file is present and generates one if not.
     *
     * @param sampleEntry The {@link SampleEntry} for which
     * @throws IOException If the initialization of the file reader or the validation/generation of the index file fails.
     */
    public static void imputeVCFFileReader(SampleEntry sampleEntry) throws IOException {
        // Check for the existence of a `.tbi.gz` index file of the input `.vcf` file.
        if (!new File(sampleEntry.vcfFile.getAbsolutePath() + ".tbi").exists()) {
            // If none is present, an index is created and written to the same directory as the input `.vcf` file.
            TabixIndex tabixIndex = IndexFactory.createTabixIndex(
                    sampleEntry.vcfFile,
                    new VCFCodec(),
                    null
            );
            tabixIndex.write(Path.of(sampleEntry.vcfFile.getAbsolutePath() + ".tbi"));
        }
        // VCFFileReader can now be initialized.
        sampleEntry.vcfFileReader = new VCFFileReader(
                sampleEntry.vcfFile,
                new File(sampleEntry.vcfFile.getAbsolutePath() + ".tbi")
        );
    }

    /**
     * Adds information about the genotype or proteoform of this sample wrt. a specific {@link FeatureEntry}.
     * <p>
     * The information is added as annotation by setting the {@link ProteoformEntry#name} or {@link AlleleEntry#name}
     * as value accessible via the key {@link FeatureEntry#name}_(PROTEOFORM|GENOTYPE).
     *
     * @param featureEntry The {@link FeatureEntry} to whose {@link AlleleEntry}/{@link ProteoformEntry} the sample should be assigned.
     * @param type         Either a {@link AlleleEntry} or {@link ProteoformEntry}.
     */
    public void associateWithType(FeatureEntry featureEntry, Object type) {
        if (type instanceof ProteoformEntry) {
            if (featureEntry.proteoforms.containsKey(((ProteoformEntry) type).name)) {
                this.annotations.put(featureEntry.name + "_PROTEOFORM", ((ProteoformEntry) type).name);
            } else {
                Logging.logWarning("Failed to assign sample " + this.name + " to proteoform "
                        + ((ProteoformEntry) type).name + " of feature "
                        + featureEntry.name + "; The specified proteoform does not exist for the feature."
                );
            }
        } else if (type instanceof AlleleEntry) {
            if (featureEntry.alleles.containsKey(((AlleleEntry) type).name)) {
                this.annotations.put(featureEntry.name + "_GENOTYPE", ((AlleleEntry) type).name);
            } else {
                Logging.logWarning("Failed to assign sample " + this.name + " to genotype "
                        + ((AlleleEntry) type).name + " of feature "
                        + featureEntry.name + "; The specified genotype does not exist for the feature."
                );
            }
        } else {
            Logging.logWarning("Failed to assign sample " + this.name + " to genotype/proteoform of feature "
                    + featureEntry.name + "; The passed type " + type.getClass().getName() + " is not supported."
            );
        }
    }

}

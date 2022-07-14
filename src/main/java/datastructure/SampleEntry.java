package datastructure;

import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;

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
     * TODO
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
     * TODO
     */
    public final HashMap<String, String> annotations = new HashMap<>();

    /**
     * Constructor of {@link SampleEntry}.
     *
     * @param vcfFileReader {@link VCFFileReader} instances initialized with the `.vcf` file of the respective sample.
     * @param name          {@link String} representing the sample name.
     */
    public SampleEntry(File vcfFile, String name) {
        this.vcfFile = vcfFile;
        this.name = name;
    }

    public static void imputeVCFFileReader(SampleEntry sampleEntry) throws IOException {
        // Check for the existence of a `.tbi.gz` index file of the input `.vcf` file.
        if (!new File(sampleEntry.vcfFile.getAbsolutePath() + ".tbi.gz").exists()) {
            // If none is present, an index is created and written to the same directory as the input `.vcf` file.
            TabixIndex tabixIndex = IndexFactory.createTabixIndex(
                    sampleEntry.vcfFile,
                    new VCFCodec(),
                    null
            );
            tabixIndex.write(Path.of(sampleEntry.vcfFile.getAbsolutePath() + ".tbi.gz"));
        }
        // VCFFileReader can now be initialized.
        sampleEntry.vcfFileReader = new VCFFileReader(
                sampleEntry.vcfFile,
                new File(sampleEntry.vcfFile.getAbsolutePath() + ".tbi.gz")
        );
    }

}

package datastructure;

import htsjdk.variant.vcf.VCFFileReader;

/**
 * Internal representation of a sample that is subject to analysis.
 * <p>
 * Pairs of instances of this class and {@link ReferenceAnalysisEntry} instances are used internally as so called 'Run
 * entries' to specify a set of analysis tasks, i.e. which sample is analyzed with respect to which specified
 * reference feature.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class SampleAnalysisEntry {

  /**
   * A {@link VCFFileReader} instance pointing to the `.vcf` file of the respective sample.
   */
  public final VCFFileReader vcfFileReader;
  /**
   * The internal name of the sample.
   */
  public final String sampleName;

  /**
   * Constructor of {@link SampleAnalysisEntry}.
   *
   * @param vcfFileReader {@link VCFFileReader} instances initialized with the `.vcf` file of the respective sample.
   * @param sampleName    {@link String} representing the sample name.
   */
  public SampleAnalysisEntry(VCFFileReader vcfFileReader, String sampleName) {
    this.vcfFileReader = vcfFileReader;
    this.sampleName = sampleName;
  }

}

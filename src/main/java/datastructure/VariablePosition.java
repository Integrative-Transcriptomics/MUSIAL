package datastructure;

import java.util.ArrayList;

/**
 * Represents and stores information about a single genomic position of one sample relative to a reference (genetic
 * sequence).
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class VariablePosition {

  /**
   * Indicates that the reference base is given, no variation in this position.
   */
  public static final char REFERENCE = '.';
  /**
   * Indicates that the reference base is given, but the position quality was low.
   */
  public static final char REFERENCE_LOW_QUAL = ':';
  /**
   * Indicates that the reference base is given, but the position coverage was low.
   */
  public static final char REFERENCE_LOW_COV = ';';
  /**
   * Represents a deletion (relative to the reference).
   */
  public static final char DELETION = '-';
  /**
   * Indicates an alternative base change to A.
   */
  public static final char ALT_A = 'A';
  /**
   * Indicates an alternative base change to A, but the quality was low.
   */
  public static final char ALT_A_LOW_QUAL = 'a';
  /**
   * Indicates an alternative base change to A, but the coverage was low.
   */
  public static final char ALT_A_LOW_COV = 'q';
  /**
   * Indicates an alternative base change to A, but the frequency was low.
   */
  public static final char ALT_A_LOW_FREQ = 'y';
  /**
   * Indicates an alternative base change to G.
   */
  public static final char ALT_G = 'G';
  /**
   * Indicates an alternative base change to G, but the quality was low.
   */
  public static final char ALT_G_LOW_QUAL = 'g';
  /**
   * Indicates an alternative base change to G, but the coverage was low.
   */
  public static final char ALT_G_LOW_COV = 'f';
  /**
   * Indicates an alternative base change to G, but the frequency was low.
   */
  public static final char ALT_G_LOW_FREQ = 'h';
  /**
   * Indicates an alternative base change to T.
   */
  public static final char ALT_T = 'T';
  /**
   * Indicates an alternative base change to T, but the quality was low.
   */
  public static final char ALT_T_LOW_QUAL = 't';
  /**
   * Indicates an alternative base change to T, but the coverage was low.
   */
  public static final char ALT_T_LOW_COV = 'r';
  /**
   * Indicates an alternative base change to T, but the frequency was low.
   */
  public static final char ALT_T_LOW_FREQ = 'z';
  /**
   * Indicates an alternative base change to C.
   */
  public static final char ALT_C = 'C';
  /**
   * Indicates an alternative base change to C, but the quality was low.
   */
  public static final char ALT_C_LOW_QUAL = 'c';
  /**
   * Indicates an alternative base change to C, but the coverage was low.
   */
  public static final char ALT_C_LOW_COV = 'v';
  /**
   * Indicates an alternative base change to C, but the frequency was low.
   */
  public static final char ALT_C_LOW_FREQ = 'd';

  /**
   * A list representing the content of the variable position, i.e. DNA bases or related characters to indicate
   * filtering values below specified thresholds (coverage, quality, frequency).
   */
  public final ArrayList<Character> content;
  /**
   * The coverage of the variable position, i.e. the number of reads spanning the position once mapped to the
   * reference.
   */
  public final double coverage;
  /**
   * The Phred scaled quality of the variable position. For more information see: <a href=https://samtools.github
   * .io/hts-specs/VCFv4.2.pdf>VCFv4.2.pdf</a>.
   */
  public final double quality;
  /**
   * A list of annotation strings of the position.
   * <p>
   *
   * @TODO: Currently not used, should be used to store SnpEff annotations once the procedure is fully implemented.
   */
  public final ArrayList<String> annotation;
  /**
   * A list storing the allele frequencies of the variable position.
   */
  public final ArrayList<Double> frequency;
  /**
   * The number of alleles of the variable position.
   */
  public int alleleCount;

  /**
   * Constructs a {@link VariablePosition} with meta-information.
   *
   * @param content    {@link Character} representing the variants base change.
   * @param coverage   The coverage of the alternative.
   * @param quality    The quality score of the alternative.
   * @param annotation Any annotations for the alternative.
   */
  public VariablePosition(char content, double coverage, double quality, String annotation,
                          double frequency) {
    this.content = new ArrayList<>();
    this.content.add(content);
    this.coverage = coverage;
    this.quality = quality;
    this.annotation = new ArrayList<>();
    this.annotation.add(annotation);
    this.frequency = new ArrayList<>();
    this.frequency.add(frequency);
    this.alleleCount = 1;
  }

  /**
   * Method to add the information about a new alternate allele to the variable position.
   *
   * @param content    {@link Character} representing the variant of the allele.
   * @param annotation {@link String} representing any annotations of the variant.
   * @param frequency  {@link Double} representing the allele frequency.
   */
  public void addAllele(char content, String annotation, double frequency) {
    this.content.add(content);
    this.annotation.add(annotation);
    this.frequency.add(frequency);
    this.alleleCount += 1;
  }

}
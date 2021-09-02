package tools;

import datastructure.VariablePositionsTable;
import exceptions.MusialIOException;
import java.io.File;
import java.util.List;
import main.Musial;
import me.tongfei.progressbar.ProgressBar;
import utility.IO;
import utility.Logging;

/**
 * Comprises a static method to generate all output files.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class OutputWriter {

  /**
   * Writes output files based on the information of a {@link VariablePositionsTable}.
   * <p>
   * Currently the following output files are generated for each reference location:
   * - Alignment of full sequences.
   * - Alignment of variant sequences.
   * - .tsv file storing the variant content for each sample and each position.
   * - .tsv file storing the variant annotations for each sample and each position.
   * - .tsv file storing the per position sum of all possible variant contents over all samples.
   * - (optional) .json file storing information for the interactive visualization extension.
   *
   * @param arguments {@link ArgumentsParser} arguments parsed from the command line.
   * @param variablePositionsTable The {@link VariablePositionsTable} of which the output should be generated,
   * @throws MusialIOException If any output file generation fails.
   */
  public static void writeOutput(ArgumentsParser arguments, VariablePositionsTable variablePositionsTable)
      throws MusialIOException {
    ProgressBar progress = Musial.buildProgress();
    try (progress) {
      List<String> referenceLocations = variablePositionsTable.getReferenceAnalysisIds();
      int isGenerateIveConfigs = arguments.isGenerateIveConfigs() ? 1 : 0;
      progress.maxHint((long) (5 + isGenerateIveConfigs) * referenceLocations.size());
      String outputDirectory = arguments.getOutputDirectory().getAbsolutePath();
      if (!outputDirectory.endsWith("/")) {
        outputDirectory += "/";
      }
      for (String referenceAnalysisId : referenceLocations) {
        IO.generateDirectory(new File(outputDirectory + "/" + referenceAnalysisId));
        // 1. Generate snv table .tsv output.
        File snvTableOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/snvTable.tsv");
        IO.generateFile(snvTableOutFile);
        IO.writeSnvTable(snvTableOutFile, referenceAnalysisId, variablePositionsTable, arguments.isAddReference());
        progress.step();
        // 2. Generate variant alignment .fasta output.
        File variantAlignmentOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/variantAlignment" +
            ".fasta");
        IO.generateFile(variantAlignmentOutFile);
        IO.writeVariantAlignment(variantAlignmentOutFile, referenceAnalysisId, variablePositionsTable, arguments.isAddReference());
        progress.step();
        // 3. Generate full alignment .fasta output.
        File fullAlignmentOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/fullAlignment.fasta");
        IO.generateFile(fullAlignmentOutFile);
        IO.writeSequenceAlignment(fullAlignmentOutFile, referenceAnalysisId, variablePositionsTable, arguments.isAddReference());
        progress.step();
        // 4. Generate snv annotations .tsv output.
        File snvAnnotationsOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/snvAnnotations.tsv");
        IO.generateFile(snvAnnotationsOutFile);
        IO.writeSnvAnnotations(snvAnnotationsOutFile, referenceAnalysisId, variablePositionsTable);
        progress.step();
        // 5. Generate per sample statistics.
        /* TODO: Not implemented yet.
        File sampleStatisticsOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/sampleStatistics" +
            ".tsv");
        IO.generateFile(sampleStatisticsOutFile);
        IO.writeSampleStatistics(fullAlignmentOutFile, referenceAnalysisId, variablePositionsTable);
        progress.step();
         */
        // 6. Generate per position statistics.
        File positionStatisticsOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/positionStatistics" +
            ".tsv");
        IO.generateFile(positionStatisticsOutFile);
        IO.writePositionStatistics(positionStatisticsOutFile, referenceAnalysisId, variablePositionsTable);
        progress.step();
        // 7. Generate IVE config files if specified by the user.
        if ( arguments.isGenerateIveConfigs() ) {
          File iveConfigOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/iveConfig.json");
          IO.generateFile(iveConfigOutFile);
          IO.writeIveConfig(iveConfigOutFile,referenceAnalysisId,variablePositionsTable,arguments);
          progress.step();
        }
      }
      progress.setExtraMessage(Logging.getDoneMessage());
    }
  }

}

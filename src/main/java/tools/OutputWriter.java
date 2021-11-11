package tools;

import datastructure.FeatureAnalysisEntry;
import datastructure.VariablePositionsTable;
import exceptions.MusialIOException;
import java.io.File;
import java.util.HashSet;
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
   * @param arguments              {@link ArgumentsParser} arguments parsed from the command line.
   * @param variablePositionsTable The {@link VariablePositionsTable} of which the output should be generated,
   * @param featureAnalysisEntries represents a genomic region of the reference data, for example a single gene, contig,
   *                               plasmid or full genome. Each such entry contains naming as well as reference sequence information.
   * @throws MusialIOException If any output file generation fails.
   */
  public static void writeOutput(ArgumentsParser arguments, VariablePositionsTable variablePositionsTable,
                                 HashSet<FeatureAnalysisEntry> featureAnalysisEntries)
      throws MusialIOException {
    ProgressBar progress = Musial.buildProgress();
    try (progress) {
      // Initialize progress counter.
      int numOutputFiles = 5;
      if (arguments.getPdInputFiles().size() != 0) {
        numOutputFiles += 3;
      }
      progress.maxHint((long) numOutputFiles * featureAnalysisEntries.size());
      // Fix output directory path.
      String outputDirectory = arguments.getOutputDirectory().getAbsolutePath();
      if (!outputDirectory.endsWith("/")) {
        outputDirectory += "/";
      }
      for (FeatureAnalysisEntry featureAnalysisEntry : featureAnalysisEntries) {
        String referenceAnalysisId = featureAnalysisEntry.identifier;
        IO.generateDirectory(new File(outputDirectory + "/" + referenceAnalysisId));
        // 1. Generate snv table .tsv output.
        File snvTableOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/snvTable.tsv");
        IO.generateFile(snvTableOutFile);
        IO.writeSnvTable(snvTableOutFile, referenceAnalysisId, variablePositionsTable);
        progress.step();
        // 2. Generate snv annotations .tsv output.
        File snvAnnotationsOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/snvAnnotations.tsv");
        IO.generateFile(snvAnnotationsOutFile);
        IO.writeSnvAnnotations(snvAnnotationsOutFile, referenceAnalysisId, variablePositionsTable);
        progress.step();
        // 3. Generate variant alignment .fasta output.
        File variantAlignmentOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/variantAlignment" +
            ".fasta");
        IO.generateFile(variantAlignmentOutFile);
        IO.writeVariantAlignment(variantAlignmentOutFile, referenceAnalysisId, variablePositionsTable);
        progress.step();
        // 4. Generate full alignment .fasta output.
        File fullAlignmentOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/fullAlignment.fasta");
        IO.generateFile(fullAlignmentOutFile);
        IO.writeFullAlignment(fullAlignmentOutFile, referenceAnalysisId, variablePositionsTable);
        progress.step();
        // 6. Generate per sample statistics.
        File sampleStatisticsOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/sampleStatistics" +
            ".tsv");
        IO.generateFile(sampleStatisticsOutFile);
        IO.writeSampleStatistics(sampleStatisticsOutFile, referenceAnalysisId, variablePositionsTable);
        progress.step();
        // 7. Generate per position statistics.
        File positionStatisticsOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/positionStatistics" +
            ".tsv");
        IO.generateFile(positionStatisticsOutFile);
        IO.writePositionStatistics(positionStatisticsOutFile, referenceAnalysisId, variablePositionsTable);
        progress.step();
        if (arguments.getPdInputFiles().containsKey(referenceAnalysisId)) {
          // 8. Generate protein data integrated results.
          File proteinDataIntegratedOutDir = new File(outputDirectory + "/" + referenceAnalysisId +
              "/proteinDataIntegratedResults");
          IO.generateDirectory(proteinDataIntegratedOutDir);
          // Generate protein alignment.
          IO.writeProteinAlignment(proteinDataIntegratedOutDir, featureAnalysisEntry, variablePositionsTable,
              arguments);
          progress.step();
          // Generate protein contact map.
          File proteinContactsMapOutFile =
              new File(proteinDataIntegratedOutDir.getAbsolutePath() + "/residueContactMap.tsv");
          IO.generateFile(proteinContactsMapOutFile);
          IO.writeProteinContactMap(proteinContactsMapOutFile, referenceAnalysisId, arguments);
          progress.step();
          // Copy .pdb file.
          File pdbFile = arguments.getPdInputFiles().get(referenceAnalysisId);
          IO.copyFile(pdbFile, new File(proteinDataIntegratedOutDir.getAbsolutePath() + "/" + pdbFile.getName()));
        } else {
          progress.stepBy(2);
        }

      }
      progress.setExtraMessage(Logging.getDoneMessage());
    }
  }

}

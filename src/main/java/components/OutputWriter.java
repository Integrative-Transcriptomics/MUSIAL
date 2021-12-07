package components;

import datastructure.FeatureAnalysisEntry;
import datastructure.VariantPositionsTable;
import exceptions.MusialBioException;
import exceptions.MusialIOException;
import java.io.File;
import java.util.HashSet;
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
   * Writes output files based on the information of a {@link VariantPositionsTable}.
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
   * @param variantPositionsTable  The {@link VariantPositionsTable} of which the output should be generated,
   * @param featureAnalysisEntries represents a genomic region of the reference data, for example a single gene, contig,
   *                               plasmid or full genome. Each such entry contains naming as well as reference sequence information.
   * @throws MusialIOException If any output file generation fails.
   */
  public static void writeOutput(ArgumentsParser arguments, VariantPositionsTable variantPositionsTable,
                                 HashSet<FeatureAnalysisEntry> featureAnalysisEntries)
      throws MusialIOException, MusialBioException {
    ProgressBar progress = Musial.buildProgress();
    try (progress) {
      // Initialize progress counter.
      int numOutputFiles = 5;
      if (arguments.getPdInputFiles().size() != 0) {
        numOutputFiles += 3;
      }
      if (arguments.isRunSnpEff()) {
        numOutputFiles += 1;
      }
      progress.maxHint((long) numOutputFiles * featureAnalysisEntries.size());
      // Fix output directory path.
      String outputDirectory = arguments.getOutputDirectory().getAbsolutePath();
      if (!outputDirectory.endsWith("/")) {
        outputDirectory += "/";
      }
      for (FeatureAnalysisEntry featureAnalysisEntry : featureAnalysisEntries) {
        String referenceAnalysisId = featureAnalysisEntry.name;
        IO.generateDirectory(new File(outputDirectory + "/" + referenceAnalysisId));
        // 1. Generate snv table .tsv output.
        File snvTableOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/variantContentTable.tsv");
        IO.generateFile(snvTableOutFile);
        IO.writeVariantsContentTable(snvTableOutFile, referenceAnalysisId, variantPositionsTable);
        progress.step();
        // 2. Generate snv annotations .tsv output.
        File snvAnnotationsOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/variantAnnotationTable" +
            ".tsv");
        IO.generateFile(snvAnnotationsOutFile);
        IO.writeVariantsAnnotationTable(snvAnnotationsOutFile, referenceAnalysisId, variantPositionsTable);
        progress.step();
        // 3. Generate sequences .fasta output.
        File sequencesOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/sequences.fasta");
        IO.generateFile(sequencesOutFile);
        IO.writeSequences(sequencesOutFile, referenceAnalysisId, variantPositionsTable);
        progress.step();
        // 4. Generate per sample statistics.
        File sampleStatisticsOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/sampleStatistics" +
            ".tsv");
        IO.generateFile(sampleStatisticsOutFile);
        IO.writeSampleStatistics(sampleStatisticsOutFile, referenceAnalysisId, variantPositionsTable);
        progress.step();
        // 5. Generate per position statistics.
        File positionStatisticsOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/positionStatistics" +
            ".tsv");
        IO.generateFile(positionStatisticsOutFile);
        IO.writePositionStatistics(positionStatisticsOutFile, referenceAnalysisId, variantPositionsTable);
        progress.step();
        if (arguments.getPdInputFiles().containsKey(referenceAnalysisId)) {
          // 6. Generate protein data integrated results.
          File proteinDataIntegratedOutDir = new File(outputDirectory + "/" + referenceAnalysisId +
              "/structureIntegratedResults");
          IO.generateDirectory(proteinDataIntegratedOutDir);
          // 6.1 Match reference nucleotide data to protein structure and allocate sample data.
          File matchingOutFile =
              new File(proteinDataIntegratedOutDir.getAbsolutePath() + "/structureAllocation.json");
          IO.generateFile(matchingOutFile);
          IO.writeStructureMatching(matchingOutFile, featureAnalysisEntry, variantPositionsTable,
              arguments);
          progress.step();
          // 6.2 Generate protein contact map.
          File proteinContactsMapOutFile =
              new File(proteinDataIntegratedOutDir.getAbsolutePath() + "/residueDistanceMap.tsv");
          IO.generateFile(proteinContactsMapOutFile);
          IO.writeProteinDistanceMap(proteinContactsMapOutFile, referenceAnalysisId, arguments);
          progress.step();
          // 6.3 Copy .pdb file.
          File pdbFile = arguments.getPdInputFiles().get(referenceAnalysisId);
          IO.copyFile(pdbFile, new File(proteinDataIntegratedOutDir.getAbsolutePath() + "/" + pdbFile.getName()));
          progress.step();
        }
        if ( arguments.isRunSnpEff() ) {
          // 7. Write SnpEff summary file.
          File snpEffSummaryOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/snpEffSummary.tsv");
          IO.generateFile(snpEffSummaryOutFile);
          IO.writeSnpEffSummary(snpEffSummaryOutFile, referenceAnalysisId, variantPositionsTable);
          progress.step();
        }
        // 8. Write run information file.
        File runInfoOutFile = new File(outputDirectory + "/" + referenceAnalysisId + "/runInfo.txt");
        IO.generateFile(runInfoOutFile);
        IO.writeRunInfo(runInfoOutFile, arguments, referenceAnalysisId);
        progress.step();
      }
      progress.setExtraMessage(Logging.getDoneMessage());
    }
  }

}

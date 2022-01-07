package components;

import datastructure.ReferenceFeatureEntry;
import datastructure.VariantContentTable;
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
   * Writes output files based on the information of a {@link VariantContentTable}.
   * <p>
   * Currently the following output files are generated for each reference feature:
   * - Alignment of full sequences.
   * - Alignment of variant sequences.
   * - .tsv file storing the variant content for each sample and each position.
   * - .tsv file storing the variant annotations for each sample and each position.
   * - .tsv file storing the per position sum of all possible variant contents over all samples.
   * - (optional) .json file storing information for the interactive visualization extension.
   *
   * @param arguments               {@link ArgumentsParser} arguments parsed from the command line.
   * @param variantContentTable     The {@link VariantContentTable} of which the output should be generated,
   * @param referenceFeatureEntries represents a genomic region of the reference data, for example a single gene,
   *                                contig,
   *                                plasmid or full genome. Each such entry contains naming as well as reference sequence information.
   * @throws MusialIOException  If any output file generation fails.
   * @throws MusialBioException If any subsequent method fails with a {@link MusialBioException}.
   */
  public static void writeOutput(ArgumentsParser arguments, VariantContentTable variantContentTable,
                                 HashSet<ReferenceFeatureEntry> referenceFeatureEntries)
      throws MusialIOException, MusialBioException {
    ProgressBar progress = Musial.buildProgress();
    // Initialize progress counter.
    int numOutputFiles = 6;
    if (arguments.getPdInputFiles().size() != 0) {
      numOutputFiles += 3;
    }
    if (arguments.isRunSnpEff()) {
      numOutputFiles += 1;
    }
    progress.maxHint((long) numOutputFiles * referenceFeatureEntries.size());
    Logging.logStatus(
        "Generating output files per reference feature (" + (long) numOutputFiles * referenceFeatureEntries.size() +
            ").");
    try (progress) {
      // Fix output directory path.
      String outputDirectory = arguments.getOutputDirectory().getAbsolutePath();
      if (!outputDirectory.endsWith("/")) {
        outputDirectory += "/";
      }
      for (ReferenceFeatureEntry referenceFeatureEntry : referenceFeatureEntries) {
        String referenceFeatureName = referenceFeatureEntry.name;
        IO.generateDirectory(new File(outputDirectory + "/" + referenceFeatureName));
        // 1. Generate snv table .tsv output.
        File snvTableOutFile = new File(outputDirectory + "/" + referenceFeatureName + "/variantContentTable.tsv");
        IO.generateFile(snvTableOutFile);
        IO.writeVariantsContentTable(snvTableOutFile, referenceFeatureName, variantContentTable);
        progress.step();

        // 2. Generate snv annotations .tsv output.
        File snvAnnotationsOutFile = new File(outputDirectory + "/" + referenceFeatureName + "/variantAnnotationTable" +
            ".tsv");
        IO.generateFile(snvAnnotationsOutFile);
        IO.writeVariantsAnnotationTable(snvAnnotationsOutFile, referenceFeatureName, variantContentTable);
        progress.step();

        // 3. Generate sequences .fasta output.
        File sequencesOutFile = new File(outputDirectory + "/" + referenceFeatureName + "/sequences.fasta");
        IO.generateFile(sequencesOutFile);
        IO.writeSequences(sequencesOutFile, referenceFeatureEntry, variantContentTable);
        progress.step();

        // 4. Generate per sample statistics.
        File sampleStatisticsOutFile = new File(outputDirectory + "/" + referenceFeatureName + "/sampleStatistics" +
            ".tsv");
        IO.generateFile(sampleStatisticsOutFile);
        IO.writeSampleStatistics(sampleStatisticsOutFile, referenceFeatureName, variantContentTable);
        progress.step();

        // 5. Generate per position statistics.
        File positionStatisticsOutFile = new File(outputDirectory + "/" + referenceFeatureName + "/positionStatistics" +
            ".tsv");
        IO.generateFile(positionStatisticsOutFile);
        IO.writePositionStatistics(positionStatisticsOutFile, referenceFeatureName, variantContentTable);
        progress.step();

        if (arguments.getPdInputFiles().containsKey(referenceFeatureName)) {
          // 6. Generate protein data integrated results.
          File proteinDataIntegratedOutDir = new File(outputDirectory + "/" + referenceFeatureName +
              "/structureIntegratedResults");
          IO.generateDirectory(proteinDataIntegratedOutDir);
          // 6.1 Match reference nucleotide data to protein structure and allocate sample data.
          File superpositionOutFile =
              new File(proteinDataIntegratedOutDir.getAbsolutePath() + "/structureSuperposition.json");
          IO.generateFile(superpositionOutFile);
          IO.writeStructureSuperposition(superpositionOutFile, referenceFeatureEntry, variantContentTable,
              arguments);
          progress.step();
          // 6.2 Generate protein contact map.
          File proteinContactsMapOutFile =
              new File(proteinDataIntegratedOutDir.getAbsolutePath() + "/residueDistanceMap.tsv");
          IO.generateFile(proteinContactsMapOutFile);
          IO.writeProteinDistanceMap(proteinContactsMapOutFile, referenceFeatureName, arguments);
          progress.step();
          // 6.3 Copy .pdb file.
          File pdbFile = arguments.getPdInputFiles().get(referenceFeatureName);
          IO.copyFile(pdbFile, new File(proteinDataIntegratedOutDir.getAbsolutePath() + "/" + pdbFile.getName()));
          progress.step();
        }

        if (arguments.isRunSnpEff()) {
          // 7. Write SnpEff summary file.
          File snpEffSummaryOutFile = new File(outputDirectory + "/" + referenceFeatureName + "/snpEffSummary.tsv");
          IO.generateFile(snpEffSummaryOutFile);
          IO.writeSnpEffSummary(snpEffSummaryOutFile, referenceFeatureName, variantContentTable);
          progress.step();
        }

        // 8. Write run information file.
        File runInfoOutFile = new File(outputDirectory + "/" + referenceFeatureName + "/runInfo.txt");
        IO.generateFile(runInfoOutFile);
        IO.writeRunInfo(runInfoOutFile, arguments, referenceFeatureName);
        progress.step();

      }
      progress.setExtraMessage(Logging.getDoneMessage());
    }
  }

}

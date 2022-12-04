package cli;

import com.google.gson.internal.LinkedTreeMap;
import components.Logging;
import components.Validation;
import exceptions.MusialException;

import java.io.File;
import java.util.ArrayList;

/**
 * Parses command line interface arguments for the `MUSIAL EXTRACT*` modules.
 * <p>
 * Used to parse and validate passed command line options or print help information to the user. Parsed arguments are
 * stored to be accessible for other components of the tool.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class ModuleExtractParameters {
    /**
     * {@link File} object specifying the input file.
     */
    public File inputFile;
    /**
     * {@link File} object specifying the output directory.
     */
    public File outputDirectory;
    /**
     * {@link ArrayList} of samples to consider.
     */
    public ArrayList<String> samples;
    /**
     * {@link ArrayList} of features to consider.
     */
    public ArrayList<String> features;
    /**
     * {@link Boolean} whether to include only SNVs.
     */
    public boolean excludeIndels;
    /**
     * {@link Boolean} whether to exclude non-variant positions.
     */
    public boolean excludeConservedPositions;
    /**
     * {@link ModuleExtractOutputModes} specifying the mode of output files generated. Either "SEQUENCE", "SEQUENCE_ALIGNED" or "TABLE".
     */
    public ModuleExtractOutputModes outputMode;
    /**
     * {@link ModuleExtractContentModes} specifying the content, i.e., either "NUCLEOTIDE" or "AMINOACID".
     */
    public ModuleExtractContentModes contentMode;
    /**
     * {@link Boolean} whether to group samples by proteoform/alleles.
     */
    public boolean grouped;


    @SuppressWarnings("FieldCanBeLocal")
    private final String EXCEPTION_PREFIX = "(Module EXTRACT Configuration)";

    /**
     * Constructor of the {@link ModuleExtractParameters} class.
     * <p>
     * Used to parse and validate command line interface input. Parsed arguments are stored - and accessible by other
     * components - via class properties.
     *
     * @param parameters {@link LinkedTreeMap} containing parameter information.
     * @throws MusialException If IO file validation fails; If CLI parameter validation fails.
     */
    public ModuleExtractParameters(LinkedTreeMap<Object, Object> parameters)
            throws MusialException {
        // Parse input file
        if (!parameters.containsKey("inputFile")) {
            throw new MusialException(EXCEPTION_PREFIX + " Missing " + Logging.colorParameter("inputFile") + "; expected path to file.");
        }
        File i = new File((String) parameters.get("inputFile"));
        if (Validation.isDirectory(new File(i.getParent()))) {
            this.inputFile = i;
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid " + Logging.colorParameter("inputFile") + " " + Logging.colorParameter((String) parameters.get("input")) + "; unable to access directory.");
        }

        // Parse output directory
        if (!parameters.containsKey("outputDirectory")) {
            throw new MusialException(EXCEPTION_PREFIX + " Missing " + Logging.colorParameter("outputDirectory") + "; expected path to directory.");
        }
        File o = new File((String) parameters.get("outputDirectory"));
        if (Validation.isDirectory(o)) {
            this.outputDirectory = o;
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid " + Logging.colorParameter("outputDirectory") + " " + Logging.colorParameter((String) parameters.get("outputDirectory")) + ".");
        }

        // Parse samples to consider
        if (parameters.containsKey("samples")) {
            //noinspection unchecked
            this.samples = (ArrayList<String>) parameters.get("samples");
        } else {
            this.samples = new ArrayList<>();
        }

        // Parse features to consider
        if (parameters.containsKey("features")) {
            //noinspection unchecked
            this.features = (ArrayList<String>) parameters.get("features");
        } else {
            this.features = new ArrayList<>();
        }

        // Parse boolean parameter whether to exclude indels.
        if (parameters.containsKey("excludeIndels")) {
            this.excludeIndels = (boolean) parameters.get("excludeIndels");
        } else {
            this.excludeIndels = false;
        }

        // Parse boolean parameter whether to include non-variant positions.
        if (parameters.containsKey("excludeConservedPositions")) {
            this.excludeConservedPositions = (boolean) parameters.get("excludeConservedPositions");
        } else {
            this.excludeConservedPositions = false;
        }

        // Parse boolean parameter whether to group samples by allele/proteoform.
        if (parameters.containsKey("grouped")) {
            this.grouped = (boolean) parameters.get("grouped");
        } else {
            this.grouped = false;
        }

        // Parse string parameter of content mode.
        if (parameters.containsKey("contentMode")) {
            try {
                this.contentMode = ModuleExtractContentModes.valueOf((String) parameters.get("contentMode"));
            } catch (IllegalArgumentException e) {
                throw new MusialException(
                        EXCEPTION_PREFIX
                                + " Unknown content mode "
                                + parameters.get("contentMode")
                );
            }
        } else {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " No content mode was specified"
            );
        }

        // Parse string parameter of output mode.
        if (parameters.containsKey("outputMode")) {
            try {
                this.outputMode = ModuleExtractOutputModes.valueOf((String) parameters.get("outputMode"));
            } catch (IllegalArgumentException e) {
                throw new MusialException(
                        EXCEPTION_PREFIX
                                + " Unknown output mode "
                                + parameters.get("outputMode")
                );
            }
        } else {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " No output mode was specified"
            );
        }
    }

}

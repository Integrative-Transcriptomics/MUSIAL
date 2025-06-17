package utility;

import com.oracle.labs.mlrg.olcut.provenance.ConfiguredObjectProvenance;
import main.Musial;
import org.tribuo.MutableDataset;
import org.tribuo.clustering.ClusterID;
import org.tribuo.clustering.ClusteringFactory;
import org.tribuo.clustering.hdbscan.HdbscanModel;
import org.tribuo.clustering.hdbscan.HdbscanTrainer;
import org.tribuo.impl.ArrayExample;
import org.tribuo.math.distance.Distance;
import org.tribuo.math.distance.L1Distance;
import org.tribuo.math.la.SGDVector;
import org.tribuo.math.neighbour.NeighboursQueryFactoryType;
import org.tribuo.math.protos.DistanceProto;
import org.tribuo.provenance.SimpleDataSourceProvenance;

import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Provides functionality for clustering data using the HDBSCAN algorithm (see <a href="https://tribuo.org/">https://tribuo.org/</a>).
 * <p>
 * The static implementation allows to re-use the trainer and model and provides methods to
 * add samples (not {@link datastructure.Sample} instances, but data points), training and
 * resetting the model, as well as retrieving clustering results.
 * <p>
 * The implementation currently only supports the clustering based on variants, i.e. each
 * feature is defined by a position and alternative base in the context of {@link datastructure.Feature.Allele}s,
 * {@link datastructure.Feature.Proteoform}s or {@link datastructure.Sample}s and considered
 * a binary feature (weight or value of 1.0). These are compared with a Manhattan distanc
 * (L1 distance).
 */
public final class Clustering {

    /**
     * Trainer for the HDBSCAN clustering algorithm with parameters:
     * <ul>
     *     <li>Minimum cluster size: 4</li>
     *     <li>Distance metric: L1 (Manhattan distance)</li>
     *     <li>Nearest neighbour query a k-d tree search</li>
     *     <li>Number of neighbours (k): 4</li>
     * </ul>
     */
    private static final HdbscanTrainer trainer = new HdbscanTrainer(4, new L1Distance(), 4, 4, NeighboursQueryFactoryType.KD_TREE);

    /**
     * The trained HDBSCAN model. Initially null until training is performed.
     */
    private static HdbscanModel model = null;

    /**
     * The clustering factory used to create clustering labels.
     */
    private static final ClusteringFactory factory = new ClusteringFactory();

    /**
     * Dataset to store the samples for clustering.
     */
    private static final MutableDataset<ClusterID> dataset
            = new MutableDataset<>(new SimpleDataSourceProvenance(Musial.runId, factory), factory);

    /**
     * List to store the names of the samples in the dataset.
     */
    private static final List<String> names = new ArrayList<>();

    /**
     * A record to represent a clustering result entry.
     *
     * @param name         The label of the sample.
     * @param label        The cluster ID assigned to the sample.
     * @param idx          The index of the sample within its cluster.
     * @param outlierScore The outlier score of the sample.
     */
    public record Entry(String name, int label, int idx, double outlierScore) {
    }

    /**
     * Resets the clustering state by clearing the dataset, labels, and model.
     */
    public static void reset() {
        dataset.clear();
        names.clear();
        model = null;
    }

    /**
     * Adds a sample to the dataset for clustering.
     *
     * @param label    The label of the sample.
     * @param variants A map of feature positions and their corresponding values.
     */
    public static void addToDataset(String label, Map<Integer, String> variants) {
        if (variants.isEmpty()) return;

        // Generate feature names by combining positions and values.
        String[] featureNames = variants.keySet().stream()
                .map(pos -> pos + variants.get(pos))
                .toArray(String[]::new);

        // Initialize feature values to 1.0.
        double[] featureValues = new double[variants.size()];
        Arrays.fill(featureValues, 1.0);

        // Add the sample to the dataset and store its label.
        dataset.add(new ArrayExample<>(factory.getUnknownOutput(), featureNames, featureValues));
        names.add(label);
    }

    /**
     * Checks if the dataset contains any samples.
     *
     * @return {@code true} if the dataset has at least one sample, {@code false} otherwise.
     */
    public static boolean hasData() {
        return dataset.size() > 0;
    }

    /**
     * Trains the HDBSCAN model using the current dataset.
     */
    public static void train() {
        Logger.getLogger(HdbscanTrainer.class.getName()).setLevel(Level.OFF);
        model = trainer.train(dataset);
    }

    /**
     * Retrieves the clustering results as a list of {@link Entry} instances.
     *
     * @return A list of clustering result entries, each containing the name, label, index, and outlier score.
     */
    public static List<Entry> getClusteringResult() {
        // Get cluster labels and outlier scores from the model.
        List<Integer> clusterLabels = model.getClusterLabels();
        List<Double> outlierScores = model.getOutlierScores();

        // Construct label re-mapping to ensure consistent cluster IDs.
        // NOTE: This is a workaround as cluster labels do not start from 1.
        Map<Integer, Integer> labelMap = new HashMap<>();
        labelMap.put(0, 0); // Outlier label

        for (int label : clusterLabels) {
            if (label != 0 && !labelMap.containsKey(label)) {
                labelMap.put(label, labelMap.size());
            }
        }

        // Map to track the size of each cluster.
        Map<Integer, Integer> clusterSizes = new HashMap<>();

        // List to store the result entries.
        List<Entry> resultEntries = new ArrayList<>();

        // Populate the result entries with clustering information.
        for (int i = 0; i < clusterLabels.size(); i++) {
            int label = labelMap.get(clusterLabels.get(i));
            clusterSizes.put(label, clusterSizes.getOrDefault(label, 0) + 1);
            resultEntries.add(new Entry(names.get(i), label, clusterSizes.get(label), outlierScores.get(i)));
        }

        return resultEntries;
    }
}

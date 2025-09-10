/**
 * This package contains classes for the execution and management of tasks within the application. It serves as the core module for
 * task-related operations.
 * <p>
 * Each implemented class should be paired with one {@link main.MusialTask} and {@link cli.CLI} implementation, and implement a
 * {@code run()} method as well as at least a {@link model.Storage} instance as field. Instances of classes implemented in the {@link op}
 * (Operations) package should be used to perform the actual atomic operations on the linked {@link model.Storage} instance.
 * <p>
 * Implemented classes are intended to be initialized in the {@link main.Musial#main} method.
 * <hr>
 * For example the {@link task.ExecutorBuild} class is responsible for executing the {@code build} task, which includes loading genomic
 * features ({@link op.FeatureLoader}), processing variant call files (VCF) ({@link op.VCFProcessor}), annotating variants
 * ({@link op.VariantAnnotator} ), inferring sequence types, and computing statistics ({@link op.StorageUpdater}), and is paired with the
 * {@link main.MusialTask#BUILD} task and {@link cli.CLIBuild} CLI implementation.
 */
package task;
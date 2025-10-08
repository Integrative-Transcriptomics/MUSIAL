/**
 * Implements atomic <b>operations</b> building the core functionality of the application.
 * <p>
 * It encapsulates the operation logic from the higher-level business logic (see {@link main}), command-line interface (see {@link cli}),
 * and underlying genomic data storage model (see {@link model}). This ensures that each operation is as self-contained as possible, can be
 * reused in different contexts, can be tested independently, and the application remains modular and maintainable.
 * <p>
 * In contrast to static utility methods (see {@link util}), operations in this package are implemented as classes that can maintain
 * intermediate processing results and states relevant for logging and modularization. Each implemented class should evolve around a linked
 * {@link model.Storage} instance and provide a public method to execute the operation.
 */
package op;
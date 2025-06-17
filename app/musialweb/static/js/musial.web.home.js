var _URL = null;
var inputReferenceSequence = null;
var inputSampleData = null;
var inputGenomicFeatures = null;
var inputFilterParameters = null;

/**
 * Class to store input reference sequence.
 */
class InputReferenceSequence {
  /**
   * Stores the reference sequence content.
   * @type {string|null}
   */
  content = null;
  /**
   * Object representing excluded positions.
   * @type {Object}
   * @todo: Entries are only removed, when users close single tags.
   */
  excludedPositions = {};
  /**
   * Object representing excluded variants at positions.
   * @type {Object}
   * @todo: Entries are only removed, when users close single tags.
   */
  excludedVariants = {};

  /**
   * Constructor of the class.
   * @constructor
   */
  constructor() {
    $("#home-input§reference-sequence").on("change", (_) =>
      this.parseContent()
    );
  }

  /**
   * Validates the exclude position value.
   *
   * @param {string} val - The exclude position value to validate.
   * @param {any} _ - Unused parameter.
   * @returns {boolean} - Returns true if the exclude position value is valid, otherwise false.
   */
  validateExcludePosition(val, _) {
    if (/^\S+\s\S+$/.test(val)) {
      let contig = val.split(" ")[0];
      let positions = val.split(" ")[1].split(",");
      return (
        /^[a-zA-Z0-9\_\-\.]+$/.test(contig) &&
        positions.every((p) => {
          return /^([0-9]+|[0-9]+\-[0-9]+)$/.test(p);
        })
      );
    } else {
      return false;
    }
  }

  /**
   * Sets `this.excludedPositions` based on the provided values.
   *
   * @todo: For some reason, one has to replace instance references ('this') with the global object reference ('inputReferenceSequence').
   *
   * @returns {void}
   */
  setExcludePosition() {
    inputReferenceSequence.excludedPositions = {};
    Metro.getPlugin("#home-input§excluded-positions", "taginput")
      .val()
      .forEach((tag) => {
        let [contig, positions] = tag.split(" ");
        if (!inputReferenceSequence.excludedPositions.hasOwnProperty(contig)) {
          inputReferenceSequence.excludedPositions[contig] = [];
        }
        for (let position of positions.split(",")) {
          if (position.includes("-")) {
            let limits = position.split("-");
            let start = parseInt(limits[0]);
            let end = parseInt(limits[1]);
            for (let i = start; i <= end; i++) {
              inputReferenceSequence.excludedPositions[contig].push(i);
            }
          } else {
            inputReferenceSequence.excludedPositions[contig].push(
              parseInt(position)
            );
          }
        }
      });
    // Make all entries unique.
    for (const [contig, entryList] of Object.entries(
      inputReferenceSequence.excludedPositions
    )) {
      inputReferenceSequence.excludedPositions[contig] = [
        ...new Set(entryList),
      ];
    }
  }

  /**
   * Validates the exclude variants value.
   *
   * @param {string} val - The exclude variants value to validate.
   * @param {any} _ - Unused parameter.
   * @returns {boolean} - Returns true if the exclude variant value is valid, otherwise false.
   */
  validateExcludeVariant(val, _) {
    if (/^\S+\s\S+$/.test(val)) {
      let contig = val.split(" ")[0];
      let entries = val.split(" ")[1].split(",");
      return (
        /^[a-zA-Z0-9\_\-\.]+$/.test(contig) &&
        entries.every((e) => {
          return /^[0-9]+\:[ACGTN]+$/.test(e);
        })
      );
    } else {
      return false;
    }
  }

  /**
   * Sets `this.excludedVariants` based on the provided values.
   *
   * @todo: For some reason, one has to replace instance references ('this') with the global object reference ('inputReferenceSequence').
   *
   * @returns {void}
   */
  setExcludeVariant() {
    inputReferenceSequence.excludedVariants = {};
    Metro.getPlugin("#home-input§excluded-variants", "taginput")
      .val()
      .forEach((tag) => {
        let [contig, entries] = tag.split(" ");
        if (!inputReferenceSequence.excludedVariants.hasOwnProperty(contig))
          inputReferenceSequence.excludedVariants[contig] = [];
        inputReferenceSequence.excludedVariants[contig].push(
          ...entries.split(",")
        );
      });
    // Make all entries unique.
    for (const [contig, entryList] of Object.entries(
      inputReferenceSequence.excludedVariants
    )) {
      inputReferenceSequence.excludedVariants[contig] = [...new Set(entryList)];
    }
  }

  /**
   * Parse the content of the reference sequence file.
   * @returns {Promise<void>} A promise that resolves when the content has been parsed.
   */
  async parseContent() {
    let referenceSequenceFile = $("#home-input§reference-sequence")[0].files[0];
    try {
      document.body.style.cursor = "wait";
      await readFile(referenceSequenceFile).then((response) => {
        this.content = response;
      });
    } catch (error) {
      displayError("Error reading reference sequence file: " + error);
    } finally {
      document.body.style.cursor = "default";
      canSubmit();
    }
  }
}

/**
 * Class for handling input sample data.
 */
class InputSampleData {
  /**
   * Stores the sample data content.
   * @type {Object}
   */
  content = {};

  /**
   * Constructor of the class.
   * @constructor
   */
  constructor() {
    $("#home-input§samples").on("change", (_) => this.parseContent());
    $("#home-input§samples-meta").on("change", (_) => this.setMetaData());
  }

  /**
   * Parses metadata for input samples and stores in `this.content`.
   */
  setMetaData() {
    try {
      document.body.style.cursor = "wait";
      $("#home-input§samples-meta").parse({
        config: {
          complete: (results, _) => {
            let header = results.data[0].slice(1);
            for (let rowData of results.data.slice(1)) {
              let sampleName = rowData[0];
              if (!this.content.hasOwnProperty(sampleName)) {
                this.content[sampleName] = {
                  vcfFile: null,
                  annotations: {},
                };
              }
              header.forEach((columnName, index) => {
                this.content[sampleName].annotations[columnName] =
                  rowData[index + 1];
              });
            }
          },
          error: (error, _) => {
            displayError("Error parsing Sample Annotation file: " + error);
          },
        },
      });
    } finally {
      document.body.style.cursor = "default";
    }
  }

  /**
   * Parses content of the input samples files and stores in `this.content`.
   * @async
   */
  async parseContent() {
    let sampleVcfFiles = $("#home-input§samples")[0].files;
    try {
      document.body.style.cursor = "wait";
      for (let sampleVcfFile of sampleVcfFiles) {
        let sampleName = sampleVcfFile.name.split(".")[0];
        if (!this.content.hasOwnProperty(sampleName)) {
          this.content[sampleName] = {
            vcfFile: null,
            annotations: {},
          };
        }
        await readFile(sampleVcfFile).then((response) => {
          this.content[sampleName].vcfFile = response;
        });
      }
    } catch (error) {
      displayError("Error reading sample vcf file: " + error);
    } finally {
      document.body.style.cursor = "default";
      canSubmit();
    }
  }

  /**
   * Retrieves all sample data content that has a VCF file attached.
   *
   * @returns {Object} An object containing the complete sample data content.
   */
  getCompleteContent() {
    return Object.fromEntries(
      Object.entries(this.content).filter(
        ([_, sampleData]) => sampleData.vcfFile !== null
      )
    );
  }
}

/**
 * Class to store input genomic features.
 */
class InputGenomicFeatures {
  /**
   *  Stores the tabulator.js object.
   */
  tabulator = null;
  /**
   * Stores the height of the tabulator object.
   */
  tabulatorHeight = 600;
  /**
   * Stores the table data.
   */
  tableData = [];
  /**
   * Stores the content of the input genomic features file.
   */
  content = null;

  /**
   * Constructor of the class.
   * @constructor
   */
  constructor() {
    this.tabulator = new Tabulator("#home-input§features-tabulator", {
      height: this.tabulatorHeight,
      minHeight: 200,
      rowHeight: 20,
      width: "auto",
      layout: "fitColumns",
      columnDefaults: {
        minWidth: 80,
      },
      selectableRows: true,
      columns: [
        {
          title: "CDS",
          field: "CDS",
          editor: (cell, onRendered, success, cancel, editorParams) => {
            cell._cell.value = !cell._cell.value;
            return false;
          },
          formatter: "toggle",
          formatterParams: {
            size: 12,
            onValue: true,
            offValue: false,
            onTruthy: true,
            onColor: "#2da47d",
            offColor: "#e4e5ed",
            clickable: true,
          },
          width: 100,
        },
        { title: "Contig", field: "Contig" },
        { title: "Source", field: "Source" },
        { title: "Type", field: "Type" },
        { title: "Start", field: "Start" },
        { title: "End", field: "End" },
        { title: "Score", field: "Score" },
        { title: "Strand", field: "Strand" },
        { title: "Frame", field: "Frame" },
        {
          title: "Attributes",
          field: "Attributes",
          width: "fitContent",
          formatter: (v) => {
            let html = ``;
            v._cell.value.split(";").forEach((attribute) => {
              let key = attribute.split("=")[0];
              let value = attribute.split("=")[1];
              html += `<code>${key}: ${value}</code> `;
            });
            return html;
          },
        },
      ],
      placeholder:
        "Choose a Reference Sequence Annotation file to fill this table.",
    });
    this.tabulator.on("rowSelectionChanged", (_) => canSubmit());
    this.tabulator.on("rowClick", (event, row) => {
      let field = event.target.getAttribute("tabulator-field");
      if (field === "CDS" || field == null) row.toggleSelect();
    });
    this.tabulator.on("rowSelected", (_) => {
      this.tabulatorHeight += 46; // TODO: This might have to be adjusted.
      this.tabulator.setHeight(this.tabulatorHeight);
      this.tabulator.redraw();
      _.freeze();
    });
    this.tabulator.on("rowDeselected", (_) => {
      this.tabulatorHeight -= 46; // TODO: This might have to be adjusted.
      this.tabulator.setHeight(this.tabulatorHeight);
      this.tabulator.redraw();
      _.unfreeze();
    });
    $("#home-input§features").on("change", (_) => this.parseContent());
    $("#home-input§features-search").on("change", (_) =>
      this.setTableFilter($("#home-input§features-search")[0].value)
    );
  }

  /**
   * Sets the table filter based on the provided value. Filter will be applied to all columns in an or-logic.
   *
   * @param {string} value
   * @returns
   */
  setTableFilter(value) {
    this.tabulator.clearFilter(true);
    if (value.trim() === "") return;
    let filters = [];
    for (let _ of this.tabulator.getColumns().map((_) => _.getField())) {
      filters.push({
        field: _,
        type: "like",
        value: value,
      });
    }
    this.tabulator.setFilter([filters]);
  }

  /**
   * Parses the content of the input genomic features file and store in `this.tableData` and `this.tabulator`.
   */
  parseContent() {
    try {
      document.body.style.cursor = "wait";
      var header = [
        "Contig",
        "Source",
        "Type",
        "Start",
        "End",
        "Score",
        "Strand",
        "Frame",
        "Attributes",
      ];
      $("#home-input§features").parse({
        config: {
          comments: "#",
          skipEmptyLines: true,
          complete: (results, _) => {
            this.content = "";
            this.tableData.length = 0;
            for (let rowData of results.data) {
              this.content += rowData.join("\t") + "\n";
              let tableEntry = { CDS: false };
              header.forEach((columnName, index) => {
                tableEntry[columnName] = rowData[index];
              });
              this.tableData.push(tableEntry);
            }
            this.tabulator.setData(this.tableData);
            canSubmit();
          },
          error: (error, _) => {
            displayError(
              "Error parsing Reference Sequence Annotation file: " + error
            );
          },
        },
      });
    } finally {
      document.body.style.cursor = "default";
    }
  }

  /**
   * Retrieves the selected features and their annotations.
   * @returns {Object} An object containing the selected features and their annotations.
   */
  getSelectedFeatures() {
    var selectedFeatures = {};
    for (let element of this.tabulator.getSelectedData()) {
      let annotations = {};
      element.Attributes.split(";").forEach((attribute) => {
        let key = attribute.split("=")[0];
        let value = attribute.split("=")[1];
        annotations[key] = value;
      });
      let featureName;
      let matchBy;
      if (annotations.hasOwnProperty("Name")) {
        featureName = annotations["Name"];
        matchBy = "match_Name";
      } else if (annotations.hasOwnProperty("ID")) {
        featureName = annotations["ID"];
        matchBy = "match_ID";
      } else if (annotations.hasOwnProperty("locus_tag")) {
        featureName = annotations["locus_tag"];
        matchBy = "match_locus_tag";
      } else {
        alert(
          "Feature " +
            element.Start +
            " ... " +
            element.End +
            " does not have a unique identifier."
        );
        continue;
      }
      selectedFeatures[featureName] = {
        [matchBy]: featureName,
        coding: element.CDS,
        annotations: annotations,
      };
    }
    return selectedFeatures;
  }
}

/**
 * Class to store filter parameters.
 */
class InputFilterParameters {
  /**
   * Stores minimal call coverage.
   * @type {number}
   */
  minimalCoverage;
  /**
   * Stores minimal call frequency.
   * @type {number}
   */
  minimalFrequency;

  /**
   * Constructor of the class.
   * @constructor
   */
  constructor() {
    // Set initial values for minimal coverage and frequency wrt. default values in input elements.
    this.minimalCoverage = 5;
    this.minimalFrequency = 90;
    // Observe changes to the call coverage input element.
    $("#home-input§call-coverage").on("change", (_) => {
      let oldValue = this.minimalCoverage;
      let newValue = undefined;
      let cancel = () => {
        displayError(
          "Minimal Variant Call Position Coverage requires positive integer; Provided " +
            newValue
        );
        $("#home-input§call-coverage")[0].value = oldValue;
      };
      try {
        newValue = parseInt($("#home-input§call-coverage")[0].value);
        if (newValue > 0) {
          this.minimalCoverage = newValue;
        } else {
          cancel();
        }
      } catch (error) {
        cancel();
      } finally {
        canSubmit();
      }
    });
    // Observe changes to the call frequency input element.
    $("#home-input§call-frequency").on("change", (_) => {
      let oldValue = this.minimalFrequency;
      let newValue = undefined;
      let cancel = () => {
        displayError(
          "Minimal Variant Call Position Frequency requires number between 0 and 100; Provided " +
            newValue
        );
        $("#home-input§call-frequency")[0].value = oldValue;
      };
      try {
        newValue = parseInt($("#home-input§call-frequency")[0].value);
        if (newValue >= 0 && newValue <= 100) {
          this.minimalFrequency = newValue;
        } else {
          cancel();
        }
      } catch (error) {
        cancel();
      } finally {
        canSubmit();
      }
    });
  }
}

/**
 * Initializes the home page.
 */
function init() {
  _URL = window.location.origin;
  inputReferenceSequence = new InputReferenceSequence();
  inputGenomicFeatures = new InputGenomicFeatures();
  inputSampleData = new InputSampleData();
  inputFilterParameters = new InputFilterParameters();
}

/**
 * Checks if the form can be submitted.
 */
function canSubmit() {
  if (
    inputReferenceSequence.content !== null &&
    Object.keys(inputGenomicFeatures.getSelectedFeatures()).length > 0 &&
    Object.keys(inputSampleData.content).length > 0 &&
    inputFilterParameters.minimalCoverage > 0 &&
    inputFilterParameters.minimalFrequency >= 0 &&
    inputFilterParameters.minimalFrequency <= 100
  ) {
    $("#home-input§submit").prop("disabled", false);
  } else {
    $("#home-input§submit").prop("disabled", true);
  }
}

/**
 * Submits the data to the server for analysis and redirects to the results page.
 */
function submit() {
  document.body.style.cursor = "wait";
  displayNotification(
    "The transmitted data is analyzed on the server. Once the analysis is complete, you will be redirected to the results page. If you close this page, the process will be aborted."
  );
  $("#home-input§submit").prop("disabled", true);
  var request = {
    minimalCoverage: inputFilterParameters.minimalCoverage,
    minimalFrequency: inputFilterParameters.minimalFrequency / 100,
    referenceSequence: inputReferenceSequence.content,
    excludedPositions: inputReferenceSequence.excludedPositions,
    excludedVariants: inputReferenceSequence.excludedVariants,
    referenceFeatures: inputGenomicFeatures.content,
    output: "",
    samples: inputSampleData.getCompleteContent(),
    features: inputGenomicFeatures.getSelectedFeatures(),
  };
  axios
    .post(_URL + "/session/start", pako.deflate(JSON.stringify(request)), {
      headers: {
        "Content-Type": "application/octet-stream",
        "Content-Encoding": "zlib",
      },
    })
    .then((response) => {
      if (responseOk(response)) {
        window.open(_URL + "/results");
      }
    })
    .catch((error) => {
      console.error(error);
      displayError(error.message);
    })
    .finally(() => {
      removeNotification();
      document.body.style.cursor = "default";
      $("#home-input§submit").prop("disabled", false);
    });
}

/**
 * Reloads the session via request to the server.
 *
 * @returns {Promise<void>} A promise that resolves when the session is reloaded.
 */
async function reloadSession() {
  displayNotification("Reloading session.");
  request = null;
  if ($("#home-input§session")[0].files.length == 0) {
    displayError(
      "You have to supply a MUSIAL session <code>.zlib</code> file in the respective input form."
    );
    removeNotification();
    return;
  }
  await readFile($("#home-input§session")[0].files[0]).then((response) => {
    request = response;
  });
  axios
    .post(_URL + "/session/reload", request, {
      headers: {
        "Content-Type": "text",
      },
    })
    .then((response) => {
      if (responseOk(response)) {
        window.open(_URL + "/results");
      }
    })
    .catch((error) => {
      console.error(error);
      displayError(error.message);
    })
    .finally(() => {
      removeNotification();
    });
}

/**
 * Initializes an example session via request to the server.
 */
function exampleSession() {
  displayNotification("Init. example session.");
  axios
    .get(_URL + "/session/example")
    .then((response) => {
      if (responseOk(response)) {
        window.open(_URL + "/results");
      }
    })
    .catch((error) => {
      console.error(error);
      displayError(error.message);
    })
    .finally(() => {
      removeNotification();
    });
}

/**
 * Read client side file to {@link String}.
 *
 * @param {String} file Client side local file path.
 * @returns {Promise<String>} File content as a {@link String}.
 */
function readFile(file) {
  return new Promise((resolve, reject) => {
    var fileReader = new FileReader();
    fileReader.onload = (event) => {
      resolve(event.target.result);
    };
    fileReader.onerror = (error) => reject(error);
    fileReader.readAsText(file);
  });
}

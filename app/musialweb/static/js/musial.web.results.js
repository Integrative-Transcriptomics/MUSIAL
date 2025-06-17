/**
 * The URL of the server to send requests to.
 */
var _URL = null;
/**
 * Global variable to store the data of the current session.
 */
var DATA = null;
/**
 * Global variable to store the records table object.
 */
var TABLE = null;
/**
 * Global variable to store the dashboard object.
 */
var DASHBOARD = null;
var DASHBOARD_OPTIONS = {};
/**
 * Global variable to indicate the active record type. One of `[samples, features, nucleotide_variants, aminoacid_variants]`.
 */
var TYPE_ACTIVE = null;
/**
 * Global variable to store the active samples. Active samples are all samples that are not filtered out by the current filters and groups.
 */
var SAMPLES_ACTIVE = [];

/**
 * Implements the records table, i.e. a Tabulator (https://tabulator.info/) object, and its functionalities.
 */
class RecordsTable {
  /**
   * Stores the tabulator.js object.
   */
  tabulator = null;
  /**
   * Stores the filters of the records table for each record type.
   */
  filter = {};
  filterSeparator = "/";
  /**
   * Stores an ECharts object to display a counts bar chart.
   */
  countsChart = {
    resizeEChart: null,
    eChartOption: countsChartOption,
    eChartPlaceholder: {
      text: "Right-click on a column header in the data record table to display counts of the corresponding column values.",
      fontSize: 11,
      fontFamily: "monospace",
      fontWeight: "bold",
      showSpinner: false,
    },
    activeColumn: null,
  };

  /**
   * Constructor of the class.
   */
  constructor() {
    this.tabulator = new Tabulator("#results-table§tabulator", {
      minHeight: "510px",
      height: "510px",
      width: "auto",
      layout: "fitColumns",
      columnDefaults: {
        minWidth: 100,
        headerTooltip: true,
        tooltip: (e, cell, onRendered) => {
          let element = document.createElement("div");
          let value = "N/A";
          if (cell.getValue() !== null)
            value = cell.getValue().replaceAll(",", " ").replaceAll(", ", " ");
          element.innerHTML =
            "<code>" + cell.getColumn().getField() + ":</code> " + value;
          return element;
        },
      },
      placeholder: "No records selected to display.",
    });
    this.tabulator.on("headerContext", (e, column) => {
      e.preventDefault();
      this._updateCountsChart(column.getField());
    });
    this.countsChart.resizeEChart = new ResizeEChart(
      "#results-table§histogram"
    );
    this.countsChart.resizeEChart.echart.setOption(
      this.countsChart.eChartOption,
      true
    );
    this.countsChart.resizeEChart.echart.showLoading(
      this.countsChart.eChartPlaceholder
    );
  }

  /**
   * Sets the counts chart data to display for the specified column. If no column is specified, the chart is cleared.
   *
   * @param {string} column The name of the column to display counts for.
   */
  _updateCountsChart = (column) => {
    if (column == undefined) {
      // If no column is specified, clear the chart and display a placeholder.
      this.countsChart.activeColumn = null;
      this.countsChart.eChartOption.title[0].text = "";
      this.countsChart.eChartOption.xAxis[0].data = [];
      this.countsChart.eChartOption.series[0].data = [];
      this.countsChart.resizeEChart.echart.setOption(
        this.countsChart.eChartOption,
        true
      );
      this.countsChart.resizeEChart.echart.showLoading(
        this.countsChart.eChartPlaceholder
      );
    } else {
      // If a column is specified, display the counts of the column values.
      let data = this._getColumnData(column, true);
      // Distinct values of the column are collected by splitting strings at ',' or ', ' symbols and adding unique values...
      let values = new Set();
      data.forEach((entry) =>
        typeof entry === "string"
          ? entry.split(/,\s|,/g).forEach((value) => values.add(value))
          : values.add(entry)
      );
      // ... and sorted in alphanumeric order.
      values = [...values].sort((v1, v2) =>
        String(v1).localeCompare(String(v2), "en", { numeric: true })
      );
      let counts = values.map(
        // The counts of the values are calculated by filtering the data for each value and counting the number of occurrences.
        (value) =>
          data.filter((entry) =>
            typeof entry === "string" ? entry.includes(value) : entry == value
          ).length
      );
      // Chart updated with the new data.
      this.countsChart.eChartOption.title[0].text =
        "Value's counts of column " + column;
      this.countsChart.eChartOption.xAxis[0].data = values;
      this.countsChart.eChartOption.series[0].data = counts;
      this.countsChart.resizeEChart.echart.setOption(
        this.countsChart.eChartOption,
        true
      );
      this.countsChart.resizeEChart.echart.hideLoading();
      this.countsChart.activeColumn = column;
    }
  };

  /**
   * Returns the data of the specified column from the table.
   *
   * @param {string} column The column identifier to retrieve data of.
   * @param {boolean} active If only records passing the current filters should be considered.
   * @returns Column data as list.
   */
  _getColumnData = (column, active) => {
    if (active)
      return this.tabulator.getData("active").map((row) => row[column]);
    else return this.tabulator.getData().map((row) => row[column]);
  };

  /**
   * Sets the records to display in the table to the specified type. The actual data is taken from the `data`-`record_<type>` entry of the cached session.
   *
   * @param {string} type One of [samples, features, nucleotide_variants, aminoacid_variants].
   */
  setRecords = (type) => {
    // Set data to according to the record type.
    let data = DATA.table_records[type];
    // Clear the table's column options.
    let columns = [];
    let columnsSelectOptions = {};
    let form_name_formatter = (cell, parameters, _) => {
      let value = "N/A";
      if (cell.getValue() !== null) {
        if (cell.getValue() == "reference")
          value = `<b style="color: #333333">R</b>`;
        else if (cell.getValue() !== undefined)
          value =
            `<b style="color: #de3c4b">` +
            cell.getValue().split(".")[0] +
            `</b>`;
        else value = cell.getValue();
      }
      return value;
    };
    // Set the columns and data of the table according to the record type.
    for (let field of data["columns"]) {
      let title = field
        .split("_")
        .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
        .join(" ");
      columns.push({
        title: title,
        field: field,
        formatter:
          field.startsWith("allele") || field.startsWith("proteoform")
            ? form_name_formatter
            : "plaintext",
      });
      columnsSelectOptions[field] = title;
    }
    columns.sort((c1, c2) => {
      if (c1.field.startsWith("name")) return -1;
      else if (c1.field.startsWith("number_of")) return -1;
      else if (c1.field.startsWith("proteoform")) return 1;
      else if (c2.field.startsWith("proteoform")) return -1;
      else if (c1.field.startsWith("allele")) return 1;
      else if (c2.field.startsWith("allele")) return -1;
      else return 0;
    });
    this.tabulator.setColumns(columns);
    this.tabulator.setData(data["rows"]);
    // Adjust Metro inputs with the new column options.
    for (let select_id of [
      "results-table§filter-field",
      "results-table§correlation-1",
      "results-table§correlation-2",
    ]) {
      Metro.getPlugin("#" + select_id, "select").data(columnsSelectOptions);
    }
    // Reset the filters, if any.
    this._refreshFilters();
    // Reset the histogram data.
    this._updateCountsChart();
  };

  /**
   * Applies the filters and groups specified in the interface to the table.
   */
  _refreshFilters = () => {
    // Remove any currently applied filters and groups.
    this.tabulator.clearFilter();
    $("#results-table§filter-list").empty();
    // Retrieve stored filters, if any.
    if (this.filter.hasOwnProperty(TYPE_ACTIVE)) {
      for (const filter of this.filter[TYPE_ACTIVE]) {
        // Create UI element per filter.
        let [field, type, value] = filter.split(this.filterSeparator);
        var item = document.createElement("div");
        item.classList.add("item");
        item.innerHTML =
          `
        <span class="label">` +
          field +
          ` ` +
          type +
          ` ` +
          value +
          `</span>
          <i class="second-action fa-duotone fa-solid fa-xmark fa-2xs"></i>`;
        $("#results-table§filter-list").append(item);
        item.childNodes[3].onclick = function () {
          TABLE.filter[TYPE_ACTIVE].delete(filter);
          item.remove();
          TABLE._refreshFilters();
        };
        // Apply filter to the table.
        this.tabulator.addFilter(field, type, value);
      }
    }
    // Update the histogram data.
    if (this.histogramActiveColumn != null)
      this._updateCountsChart(this.histogramActiveColumn);
    // Update the active samples.
    if (TYPE_ACTIVE == "samples") {
      let _ = this._getColumnData("name", true);
      if (_.length == DATA.table_records["samples"].rows.length) _ = [];
      let difference = _.filter((x) => !SAMPLES_ACTIVE.includes(x)).concat(
        SAMPLES_ACTIVE.filter((x) => !_.includes(x))
      );
      if (difference.length > 0) {
        setSamplesActive(_, () => {
          updateDashboardOption("features");
          updateDashboardOption("nucleotide_variants");
          updateDashboardOption("samples");
        });
        DASHBOARD.setContent("samples");
      }
    }
  };

  /**
   * Adds a filter tag in the correct format to the interface element.
   */
  addFilter = () => {
    if (this.filter[TYPE_ACTIVE] == undefined)
      this.filter[TYPE_ACTIVE] = new Set();
    this.filter[TYPE_ACTIVE].add(
      [
        $("#results-table§filter-field")[0].value,
        $("#results-table§filter-type")[0].value,
        $("#results-table§filter-value")[0].value,
      ].join(this.filterSeparator)
    );
    this._refreshFilters();
  };

  /**
   * Downloads the table data in CSV format.
   */
  saveCsv = () => {
    this.tabulator.download(
      "csv",
      TYPE_ACTIVE + "_table_" + Date.now() + ".csv"
    );
  };
}

/**
 * Class to store and construct EChart instances with an observer for resizing.
 */
class ResizeEChart {
  /**
   * Stores the EChart instance.
   */
  echart = null;
  /**
   * Stores the ResizeObserver instance.
   */
  _observer = null;

  /**
   * Constructor of the class.
   *
   * @param {string} div The id of the div element to create the EChart instance in.
   */
  constructor(div) {
    let chartDiv = $(div)[0];
    this.echart = echarts.init(chartDiv, {
      devicePixelRatio: 2,
      renderer: "canvas",
      width: "auto",
      height: "auto",
    });
    this._observer = new ResizeObserver((entries) => {
      this.echart.resize({
        width: entries[0].width,
        height: entries[0].height,
      });
    });
    this._observer.observe(chartDiv);
  }
}

/**
 * Handles the visualizations of data in the current session depending on the type of records displayed.
 */
class Dashboard {
  /**
   * Instance of the `Plot` class to store the EChart instance.
   */
  plot = null;

  /**
   * Constructor of the class.
   *
   * @param {string} div The id of the div element, passed to the `Plot` class.
   */
  constructor(div) {
    this.plot = new ResizeEChart(div);
  }

  /**
   * Sets the content of the dashboard depending on the type of records displayed.
   *
   * @param {string} type Record type; one of `[samples, features, nucleotide_variants, aminoacid_variants]`.
   */
  setContent = (type) => {
    switch (type) {
      case "samples":
        this._setSamplesContent();
        break;
      case "features":
        this._setFeaturesContent();
        break;
      case "nucleotide_variants":
        this._setNVariantContent();
        break;
      case "aminoacid_variants":
        this._setAVariantContent();
        break;
      default:
        break;
    }
  };

  /**
   * Sets the content of the dashboard to display the samples data.
   *
   * This function is separated in order to allow tracking of events on the dashboard, specific to the samples data.
   */
  _setSamplesContent = () => {
    this.plot.echart.off("click");
    this.plot.echart.setOption(DASHBOARD_OPTIONS["samples"], true);
  };

  /**
   * Sets the content of the dashboard to display the features data.
   *
   * This function is separated in order to allow tracking of events on the dashboard, specific to the features data.
   */
  _setFeaturesContent = () => {
    this.plot.echart.off("click");
    this.plot.echart.setOption(DASHBOARD_OPTIONS["features"], true);
  };

  /**
   * Sets the content of the dashboard to display the nucleotide variants data.
   *
   * This function is separated in order to allow tracking of events on the dashboard, specific to the nucleotide variants data.
   */
  _setNVariantContent = () => {
    // Remove any previous event listeners.
    this.plot.echart.off("click");
    // Add a new event listener: On click on a feature label, dispatch a data zoom action to zoom into the feature's region.
    this.plot.echart.on("click", (params) => {
      if (params.seriesId.startsWith("dashboard.NVariants.featureTrack."))
        DASHBOARD.plot.echart.dispatchAction({
          type: "dataZoom",
          dataZoomIndex: 1,
          startValue: params.data._zoomStart,
          endValue: params.data._zoomEnd,
        });
    });
    this.plot.echart.setOption(DASHBOARD_OPTIONS["nucleotide_variants"], true);
  };

  /**
   * Sets the content of the dashboard to display the aminoacid variants data.
   *
   * This function is separated in order to allow tracking of events on the dashboard, specific to the aminoacid variants data.
   */
  _setAVariantContent = () => {
    this.plot.echart.off("click");
    this.plot.echart.setOption(DASHBOARD_OPTIONS["aminoacid_variants"], true);
  };
}

/**
 * Handles the parameter selection and download of sequences in fasta format.
 */
class SequenceDownloadPopup {
  /**
   * Configuration object for the Swal2 popup that is displayed to users for parameter selection.
   */
  swalConfiguration = {
    customClass: {
      title: "swal2-style",
      confirmButton: "swal2-style",
      cancelButton: "swal2-style",
      htmlContainer: "swal2-style",
    },
    html:
      `
    <h5 class="text-upper">download sequences in fasta format - please specify the following parameters and proceed with download.</h5>
    <br>
      <div class="grid">
        <div class="row text-left">
          <div class="cell-4 m-2">
            <p class="text-upper">gene to download sequences for</p>
          </div>
          <div class="cell-4 m-2">
            <select
              id="feature-select"
              class="input-small"
              data-role="select"
              data-filter="false"
            >` +
      DATA.table_records["features"].rows
        .map(
          (row) => '<option value="' + row.name + '">' + row.name + "</option>"
        )
        .join("") +
      `</select>
          </div>
        </div>

        <div class="row text-left">
          <div class="cell-4 m-2">
            <p class="text-upper">sequence content</p>
          </div>
          <div class="cell-4 m-2">
            <select
              id="content-select"
              class="input-small"
              data-role="select"
              data-filter="false"
            >
              <option value="nucleotide">Nucleotide</option>
              <option value="aminoacid">Aminoacid</option>
            </select>
          </div>
        </div>

        <div class="row text-left">
          <div class="cell-4 m-2">
          <p class="text-upper">sequence compilation options</p>
          </div>
          <div class="cell-2 m-2">
            <input
              id="align-select"
              type="checkbox"
              checked
              data-role="checkbox"
              data-caption="ALIGN SEQUENCES"
              data-material="true"
            />
          </div>
          <div class="cell-2 m-2">
            <input
              id="merge-select"
              type="checkbox"
              checked
              data-role="checkbox"
              data-caption="MERGE SAMPLES"
              data-material="true"
            />
          </div>
          <div class="cell-2 m-2">
            <input
              id="conserved-select"
              type="checkbox"
              checked
              data-role="checkbox"
              data-caption="INCLUDE CONSERVED POSITIONS"
              data-material="true"
            />
          </div>
        </div>

        <div class="row text-left">
          <div class="cell-4 m-2">
            <p class="text-upper"><small>include or exclude listed samples (default/empty is all samples)</small></p>
          </div>
          <div class="cell-2 m-2">
            <button data-role="hint" data-hint-text="COPY FILTERED SAMPLES FROM TABLE" data-hint-position="right" class="ml-1 button small rounded" onclick="Metro.getPlugin('#samples-select','taginput').val(SAMPLES_ACTIVE)"><i class="fa-duotone fa-paste"></i></button>
          </div>
          <div class="cell-4 m-2">
            <input
              id="samples-mode-select"
              type="checkbox"
              checked
              data-role="checkbox"
              data-caption="INCLUDE (CHECKED) OR EXCLUDE (UNCHECKED) LISTED SAMPLES"
              data-material="true"
            />
          </div>
        </div>

        <div class="row text-left">
          <div class="cell-10 m-2">
            <input id="samples-select" class="input-small" type="text" data-role="taginput" style="overflow-y: hide"/>
          </div>
        </div>
      </div>`,
    width: "100vw",
    padding: "5em",
    position: "top",
    showCancelButton: true,
    grow: true,
    heightAuto: false,
    cancelButtonColor: "#fe4848cc",
    cancelButtonText: "Cancel",
    confirmButtonColor: "#39c093cc",
    confirmButtonText: "Download",
    color: "#747474",
    background: "transparent",
    backdrop: `
      #fafafcee
      left top
      no-repeat
    `,
  };

  /**
   * Constructor of the class.
   */
  constructor() {
    Swal.fire(this.swalConfiguration).then((result) => {
      if (result.isConfirmed) {
        var request = {
          feature: $("#feature-select")[0].value,
          content: $("#content-select")[0].value,
          align: $("#align-select").is(":checked"),
          merge: $("#merge-select").is(":checked"),
          conserved: $("#conserved-select").is(":checked"),
          samples: [],
        };
        let include = $("#samples-mode-select").is(":checked");
        let selectedSamples = Metro.getPlugin(
          "#samples-select",
          "taginput"
        ).val();
        let allSamples = DATA.table_records["samples"].rows.map(
          (row) => row.name
        );
        if (selectedSamples.length == 0) {
          request["samples"] = allSamples;
        } else if (include) {
          request["samples"] = selectedSamples;
        } else {
          request["samples"] = allSamples.filter((sampleName) => {
            return !selectedSamples.includes(sampleName);
          });
        }
        displayNotification(
          "Request sent. Your sequence data will be downloaded automatically as soon as processing is complete. If you close this page, the process will be canceled."
        );
        axios
          .post(
            _URL + "/download/sequences",
            pako.deflate(JSON.stringify(request)),
            {
              headers: {
                "Content-Type": "application/octet-stream",
                "Content-Encoding": "zlib",
              },
              responseType: "blob",
            }
          )
          .then((response) => {
            if (responseOk(response)) {
              downloadBlob(
                response.data,
                request.feature + "_sequences_" + Date.now() + ".fasta"
              );
            }
          })
          .catch((error) => {
            console.error(error);
            displayError(error.message + " " + error.response.data);
          })
          .finally(() => {
            removeNotification();
          });
      }
    });
  }
}

/**
 * Handles the parameter slection and download of variants in tsv format.
 */
class VariantsDownloadPopup {
  /**
   * Configuration object for the Swal2 popup that is displayed to users for parameter selection.
   */
  swalConfiguration = {
    customClass: {
      title: "swal2-style",
      confirmButton: "swal2-style",
      cancelButton: "swal2-style",
      htmlContainer: "swal2-style",
    },
    html:
      `
    <h5 class="text-upper">download variants in tsv format - please specify the following parameters and proceed with download.</h5>
    <br>
      <div class="grid">
        <div class="row text-left">
          <div class="cell-4 m-2">
            <p class="text-upper">gene to download variants for</p>
          </div>
          <div class="cell-4 m-2">
            <select
              id="feature-select"
              class="input-small"
              data-role="select"
              data-filter="false"
            >` +
      DATA.table_records["features"].rows
        .map(
          (row) => '<option value="' + row.name + '">' + row.name + "</option>"
        )
        .join("") +
      `</select>
          </div>
        </div>

        <div class="row text-left">
          <div class="cell-4 m-2">
            <p class="text-upper">variants content</p>
          </div>
          <div class="cell-4 m-2">
            <select
              id="content-select"
              class="input-small"
              data-role="select"
              data-filter="false"
            >
              <option value="nucleotide">Nucleotide</option>
              <option value="aminoacid">Aminoacid</option>
            </select>
          </div>
        </div>

        <div class="row text-left">
          <div class="cell-4 m-2">
            <p class="text-upper"><small>include or exclude listed samples (default/empty is all samples)</small></p>
          </div>
          <div class="cell-2 m-2">
            <button data-role="hint" data-hint-text="COPY FILTERED SAMPLES FROM TABLE" data-hint-position="right" class="ml-1 button small rounded" onclick="Metro.getPlugin('#samples-select','taginput').val(SAMPLES_ACTIVE)"><i class="fa-duotone fa-paste"></i></button>
          </div>
          <div class="cell-4 m-2">
            <input
              id="samples-mode-select"
              type="checkbox"
              checked
              data-role="checkbox"
              data-caption="INCLUDE (CHECKED) OR EXCLUDE (UNCHECKED) LISTED SAMPLES"
              data-material="true"
            />
          </div>
        </div>

        <div class="row text-left">
          <div class="cell-10 m-2">
            <input id="samples-select" class="input-small" type="text" data-role="taginput" style="overflow-y: hide"/>
          </div>
        </div>
      </div>`,
    width: "100vw",
    padding: "5em",
    position: "top",
    showCancelButton: true,
    grow: true,
    heightAuto: false,
    cancelButtonColor: "#fe4848cc",
    cancelButtonText: "Cancel",
    confirmButtonColor: "#39c093cc",
    confirmButtonText: "Download",
    color: "#747474",
    background: "transparent",
    backdrop: `
      #fafafcee
      left top
      no-repeat
    `,
  };

  /**
   * Constructor of the class.
   */
  constructor() {
    Swal.fire(this.swalConfiguration).then((result) => {
      if (result.isConfirmed) {
        var request = {
          feature: $("#feature-select")[0].value,
          content: $("#content-select")[0].value,
          samples: [],
        };
        let include = $("#samples-mode-select").is(":checked");
        let selectedSamples = Metro.getPlugin(
          "#samples-select",
          "taginput"
        ).val();
        let allSamples = DATA.table_records["samples"].rows.map(
          (row) => row.name
        );
        if (selectedSamples.length == 0) {
          request["samples"] = allSamples;
        } else if (include) {
          request["samples"] = selectedSamples;
        } else {
          request["samples"] = allSamples.filter((sampleName) => {
            return !selectedSamples.includes(sampleName);
          });
        }
        displayNotification(
          "Request sent. Your variants table data will be downloaded automatically as soon as processing is complete. If you close this page, the process will be canceled."
        );
        axios
          .post(
            _URL + "/download/variants",
            pako.deflate(JSON.stringify(request)),
            {
              headers: {
                "Content-Type": "application/octet-stream",
                "Content-Encoding": "zlib",
              },
              responseType: "blob",
            }
          )
          .then((response) => {
            if (responseOk(response)) {
              downloadBlob(
                response.data,
                request.feature + "_variants_" + Date.now() + ".tsv"
              );
            }
          })
          .catch((error) => {
            console.error(error);
            displayError(error.message + " " + error.response.data);
          })
          .finally(() => {
            removeNotification();
          });
      }
    });
  }
}

/**
 * Initializes the results page.
 */
async function init() {
  // Initialize global variables.
  _URL = window.location.origin;
  TABLE = new RecordsTable();
  DASHBOARD = new Dashboard("#results-charts§echart");
  // Request session data from the server.
  axios
    .get(_URL + "/session/data")
    .then((response) => {
      if (responseOk(response)) {
        DATA = response.data;
        console.log(DATA);
        // Update record type dependent interface elements.
        /*
         * Samples:
         * - Update the UI inputs with the samples records column names for color selection.
         */
        $("#results-chart§control-samples§color")[0].innerHTML = "";
        for (const category of [
          "table_filter",
          ...DATA.table_records.samples.columns,
        ]) {
          if (category != "name") {
            var _ = document.createElement("option");
            _.value = category;
            _.text = category.split("_").join(" ");
            $("#results-chart§control-samples§color")[0].append(_);
          }
        }
        /*
         * Features:
         * - Update the UI inputs with the features records row names for feature selection
         */
        $("#results-chart§control-features§feature")[0].innerHTML = "";
        DATA.table_records.features.rows.forEach((row) => {
          var _ = document.createElement("option");
          _.value = row.name;
          _.text = row.name;
          $("#results-chart§control-features§feature")[0].append(_);
        });
        // For each available record type...
        for (let record_type of Object.keys(DATA.table_records)) {
          // Display the number of rows in the UI.
          $("#results-control§show-" + record_type).html(
            record_type.replace("_", " ") +
              ` <span class="badge">` +
              DATA.table_records[record_type].rows.length +
              `</span>`
          );
          // Generate the (initial) dashboard options for each record type. @Todo: Exclude this step from the loop and add all available dashboard ids.
          updateDashboardOption(record_type);
        }

        // Display the samples records by default.
        displayRecords("samples");
        DASHBOARD.setContent("samples");
      }
    })
    .catch((error) => {
      console.error(error);
      displayError(error.message + " " + error.response.data);
    });
}

/**
 * Transfers user interaction on the interface to the `RecordsTable` instance.
 *
 * This method does not update any data.
 *
 * @param {string} type One of `samples`, `features`, `nucleotide_variants` or `aminoacid_variants`.
 */
async function displayRecords(type) {
  // TABLE.saveFilters(); TODO
  TYPE_ACTIVE = type;
  // Indicate the active record type.
  $(".results-control-button").removeClass("selected-category");
  $("#results-control§show-" + TYPE_ACTIVE).addClass("selected-category");
  // Show the correct chart control interface.
  $(".chartcontrol").hide();
  $("#results-chart§control-" + TYPE_ACTIVE).show();
  // Set the table records and dashboard option of the active type.
  await TABLE.setRecords(TYPE_ACTIVE);
  await DASHBOARD.setContent(TYPE_ACTIVE);
}

/**
 * Function to set the active samples for the current session. This interface allows to execute a callback function after the samples have been set.
 *
 * @param {Array} value The list of samples to set as active.
 * @param {Function} callback The function to call after the samples have been set.
 */
function setSamplesActive(value, callback) {
  SAMPLES_ACTIVE = value;
  if (callback !== undefined) callback();
}

/**
 * Updates the ECharts option object of the dashboard of the specified ID.
 *
 * The actual objects are generated by functions implemented in `musial.web.dashboard.js`. Current options for ID are:
 * - `samples` Per sample per feature variant count heatmap and tSNE embedding of sample to sample distances.
 * - `features` Per feature per sample variant count heatmap and allele network of selected feature alleles.
 * - `nucleotide_variants` Per variable position view of nucleotide variants.
 * - `aminoacid_variants` Not implemented.
 *
 * @param {string} id See list of available IDs above.
 */
async function updateDashboardOption(id) {
  switch (id) {
    case "samples":
      DASHBOARD_OPTIONS[id] = getDashboardOptionSamples();
      return;
    case "features":
      DASHBOARD_OPTIONS[id] = getDashboardOptionFeatures();
      return;
    case "nucleotide_variants":
      DASHBOARD_OPTIONS[id] = getDashboardOptionNVariants();
      return;
    case "aminoacid_variants":
      DASHBOARD_OPTIONS[id] = getDashboardOptionAVariants();
      return;
    default:
      return;
  }
}

/**
 * Re-computes the EChart option for the dashboard of the specified ID and sets it as the current dashboard content.
 *
 * The ID must be a key of the `DASHBOARD_OPTIONS` object and has to be initialized in the `init` function of this script. Available values are documented in the `updateDashboardOption` function.
 *
 * @param {string} id The ID of the dashboard option to recompute.
 */
async function refreshDashboard(id) {
  if (DASHBOARD_OPTIONS.hasOwnProperty(id)) {
    updateDashboardOption(id);
    DASHBOARD.setContent(id);
  }
}

/**
 * Transfers user interaction on the interface to the server for the download of the current user session data.
 */
function sendRequestDownloadSession() {
  displayNotification(
    "Request has been sent to the server. Your session data will be downloaded automatically as soon as processing is complete. Please do not close this page."
  );
  axios
    .get(_URL + "/session/save")
    .then((response) => {
      if (responseOk(response)) {
        downloadBlob(
          response.data,
          "MUSIALweb_Session_" + Date.now() + ".zlib"
        );
      }
    })
    .catch((error) => {
      console.error(error);
      displayError(error.message + " " + error.response.data);
    })
    .finally(() => {
      removeNotification();
    });
}

/**
 * TODO
 */
function sendRequestDownloadSequenceTyping() {
  displayNotification(
    "Request has been sent to the server. The sequence typing will be downloaded automatically as soon as processing is complete. Please do not close this page."
  );
  axios
    .get(_URL + "/download/sequence_typing")
    .then((response) => {
      if (responseOk(response)) {
        downloadBlob(response.data, "sequence_typing_" + Date.now() + ".tsv");
      }
    })
    .catch((error) => {
      console.error(error);
      displayError(error.message + " " + error.response.data);
    })
    .finally(() => {
      removeNotification();
    });
}

/**
 * Transfers user interaction on the interface to the server in order to compute the correlation between two columns of the records table.
 */
function sendRequestComputeCorrelation() {
  var request = {
    col1: Metro.getPlugin("#results-table§correlation-1", "select").val(),
    col2: Metro.getPlugin("#results-table§correlation-2", "select").val(),
    test_type: Metro.getPlugin(
      "#results-table§correlation-test",
      "select"
    ).val(),
    entry_type: TYPE_ACTIVE,
  };
  axios
    .post(
      _URL + "/compute/correlation",
      pako.deflate(JSON.stringify(request)),
      {
        headers: {
          "Content-Type": "application/octet-stream",
          "Content-Encoding": "zlib",
        },
      }
    )
    .then((response) => {
      if (responseOk(response)) {
        $("#results-table§correlation-result").html(
          `Test Stat. Value: ` +
            response.data.t +
            `&nbsp;&nbsp;&nbsp;P Value: ` +
            response.data.p
        );
      }
    })
    .catch((error) => {
      console.error(error);
      displayError(error.message + " " + error.response.data);
    });
}

/**
 * Global ECharts option specifications to display text.
 */
const textOption = {
  textStyle: {
    color: "#333333",
    fontWeight: "lighter",
    fontFamily: "monospace",
    fontSize: 11,
  },
};
/**
 * Global ECharts option specifications to display titles.
 */
const titleOption = {
  left: "center",
  textStyle: {
    fontFamily: "monospace",
    fontSize: 13,
  },
};
/**
 * Global ECharts option specifications to display axis.
 */
const axisOption = {
  axisLabel: {
    fontFamily: "monospace",
    fontSize: 13,
    fontWeight: "lighter",
  },
  axisTick: {
    alignWithLabel: true,
  },
  nameLocation: "middle",
  nameTextStyle: { fontWeight: "lighter", fontSize: 12 },
  nameGap: 30,
};
/**
 * Global ECharts option specifications to display labels.
 */
const labelOption = {
  label: {
    show: true,
    fontFamily: "monospace",
    fontSize: 11,
    color: "#fafafc",
  },
};
/**
 * Global ECharts option specifications to display tooltips.
 */
const tooltipOption = {
  show: true,
  trigger: "item",
  backgroundColor: "#fafafc",
  borderColor: "transparent",
  ...textOption,
};
/**
 * Global ECharts option specifications to display zoom sliders.
 */
const sliderOption = {
  showDetail: false,
  showDataShadow: false,
  backgroundColor: "transparent",
  fillerColor: "rgba(155, 166, 189, 0.2)",
  borderColor: "#cbd0e0",
  moveHandleStyle: {
    color: "#cbd0e0",
  },
  emphasis: {
    moveHandleStyle: {
      color: "#6d81ad",
    },
  },
};
/**
 * ECharts option specification to display a histogram of the records table.
 */
const countsChartOption = {
  grid: [
    {
      id: "RecordsTable.Histogram.Grid",
      show: false,
      top: 40,
      bottom: 20,
      left: 40,
      right: 20,
      containLabel: true,
    },
  ],
  title: [
    {
      id: "RecordsTable.Histogram.Title",
      top: 0,
      left: "center",
      text: "",
      ...titleOption,
    },
  ],
  xAxis: [
    {
      id: "RecordsTable.Histogram.XAxis",
      type: "category",
      data: [],
      name: "Value",
      ...axisOption,
      gridIndex: 0,
    },
  ],
  yAxis: [
    {
      id: "RecordsTable.Histogram.YAxis",
      type: "value",
      name: "Count",
      ...axisOption,
      gridIndex: 0,
    },
  ],
  tooltip: {
    ...tooltipOption,
    formatter: (params, ticket, callback) => {
      return (
        "Value <code>" +
        params.name +
        "</code> occurrs <code>" +
        params.value +
        "</code> times."
      );
    },
  },
  series: [
    {
      id: "RecordsTable.Histogram.Series",
      type: "bar",
      data: [],
      xAxisIndex: 0,
      yAxisIndex: 0,
      barWidth: "60%",
      itemStyle: {
        color: "#333333",
      },
    },
  ],
};
const colorMapOrRd = ["#fee8c8", "#fdbb84", "#e34a33"].reverse();

/**
 * Constructs the ECharts option specification used to fill the dashboard for the samples record type.
 *
 * @returns {Object} ECharts option specification to display via the `Dashboard` class defined in the `musial.web.results.js` module.
 */
function getDashboardOptionSamples() {
  var _feature_names = Object.keys(DATA.storage.features);
  var _sample_names_filtered =
    SAMPLES_ACTIVE.length > 0
      ? SAMPLES_ACTIVE
      : Object.keys(DATA.storage.samples);
  var visualMap = []; // Stores the visual map specifications for the ECharts option.
  var legend = []; // Stores the legend specifications for the ECharts option.

  // Construct series for samples distance manifold.
  // The coloring state selected by the user.
  var _manifold_color_by = $("#results-chart§control-samples§color")[0].value;
  // The type of coloring state.
  var _manifold_color_type = _manifold_color_by.startsWith("number_of_")
    ? "continuous"
    : "categorical";
  // Object used to assign colors to samples based on the selected coloring state.
  var _manifold_colors = {
    get: (value) => {
      if (!_manifold_colors._map.hasOwnProperty(value)) {
        if (value == "Pass") _manifold_colors._map[value] = "#de3c4b";
        else if (value == "Other") _manifold_colors._map[value] = "#333333";
        else {
          _manifold_colors._map[value] =
            _manifold_colors._colors[_manifold_colors._c];
          _manifold_colors._c += 1;
          if (_manifold_colors._c >= _manifold_colors._colors.length)
            _manifold_colors._c = 0;
        }
      }
      return _manifold_colors._map[value];
    },
    _map: {},
    _colors: [
      "#333333",
      "#a6cee3",
      "#1f78b4",
      "#b2df8a",
      "#33a02c",
      "#fb9a99",
      "#e31a1c",
      "#fdbf6f",
      "#ff7f00",
      "#cab2d6",
      "#6a3d9a",
    ],
    _c: 0,
  };
  // Noise to add to the tSNE projection to avoid overlapping points.
  var _manifold_noise = () =>
    DATA.meta_components.sample_distance_manifold.mean_nnd *
    (Math.random() - 0.5) *
    0.1;
  var _manifold_values = [];
  var _manifold_series_idx = 0;
  // Series specification for the samples distance manifold.
  var _manifold_series = {};
  Object.keys(DATA.meta_components.sample_distance_manifold.samples).forEach(
    (sample_name) => {
      let series_id = "N/A";
      if (_manifold_color_by == "" || _manifold_color_by == "table_filter") {
        series_id = SAMPLES_ACTIVE.includes(sample_name) ? "Pass" : "Other";
      } else if (_manifold_color_by.startsWith("allele")) {
        series_id =
          DATA.storage.samples[sample_name].allele[_manifold_color_by.slice(7)];
      } else if (
        _manifold_color_by.startsWith("proteoform") &&
        DATA.storage.features[_manifold_color_by.slice(11)].type == "coding"
      ) {
        series_id =
          DATA.storage.samples[sample_name].proteoform[
            _manifold_color_by.slice(11)
          ];
      } else if (_manifold_color_type == "categorical") {
        series_id = DATA.storage.samples[sample_name].info[_manifold_color_by];
      } else if (_manifold_color_type == "continuous") {
        _manifold_values.push(
          parseInt(DATA.storage.samples[sample_name].info[_manifold_color_by])
        );
      }
      let series_color = _manifold_colors.get(series_id);
      let projection =
        DATA.meta_components.sample_distance_manifold.samples[sample_name]
          .projection;
      if (!_manifold_series.hasOwnProperty(series_id)) {
        _manifold_series[series_id] = {
          id: "dashboard.Samples.manifold." + series_id,
          name: series_id,
          type: "scatter",
          xAxisIndex: 2,
          yAxisIndex: 2,
          data: [],
          z: 0,
          tooltip: {
            formatter: (params, ticket, callback) => {
              return _formatSampleTooltip(params.data.name);
            },
          },
          symbolSize: 5,
          itemStyle: {},
        };
        _manifold_series_idx += 1;
        if (_manifold_color_type == "categorical")
          _manifold_series[series_id].itemStyle.color = series_color;
      }
      _manifold_series[series_id].data.push({
        name: sample_name,
        value: [
          projection[0] + _manifold_noise(),
          projection[1] + _manifold_noise(),
          _manifold_color_type == "continuous"
            ? parseFloat(
                DATA.storage.samples[sample_name].info[_manifold_color_by]
              )
            : NaN,
        ],
      });
    }
  );
  if (_manifold_color_type == "continuous") {
    _manifold_values.sort((a, b) => a - b);
    var _manifold_lim =
      _manifold_values
        .slice(0, Math.round(_manifold_values.length * 0.95))
        .pop() + 1;
    visualMap.push({
      type: "continuous",
      seriesIndex: [...Array(_manifold_series_idx).keys()],
      min: 0,
      max: _manifold_lim,
      range: [0, _manifold_lim],
      color: colorMapOrRd,
      orient: "horizontal",
      top: 25,
      left: "56%",
      itemWidth: 10,
      text: ["> " + (_manifold_lim - 1), "Count: 0"],
      ...textOption,
    });
  } else if (_manifold_color_type == "categorical") {
    legend.push({
      show: true,
      type: "scroll",
      top: 22,
      left: "53%",
      width: "40%",
      data: Object.keys(_manifold_series).map((_) => {
        return {
          name: _,
        };
      }),
      itemGap: 25,
      ...textOption,
    });
  }

  // Construct series for variants counts.
  var _counts_series = [
    {
      id: "dashboard.Samples.counts",
      type: "heatmap",
      xAxisIndex: 0,
      yAxisIndex: 0,
      data: DATA.meta_components.variant_counts_clustering.counts
        .map((_) => [
          _[1], // Feature
          _[0], // Sample
          _[2], // Count
        ])
        .filter((_) => _[2] > 0),
      z: 0,
      tooltip: {
        formatter: (params, ticket, callback) => {
          let _ = _formatSampleTooltip(params.data[1], params.data[0]);
          _ += _formatFeatureTooltip(params.data[0]);
          return _;
        },
        position: "right",
        showDelay: 60,
        hideDelay: 20,
      },
      progressive: 2000,
      progressiveThreshold: 1000,
    },
    {
      id: "dashboard.Samples.marginal",
      type: "bar",
      xAxisIndex: 1,
      yAxisIndex: 1,
      data: Object.entries(
        DATA.meta_components.variant_counts_clustering.counts.reduce((A, c) => {
          if (A.hasOwnProperty(c[0])) A[c[0]] += c[2];
          else A[c[0]] = c[2];
          return A;
        }, {})
      ).map((_) => [_[1], _[0]]),
      itemStyle: {
        color: "#333333",
      },
      tooltip: {
        formatter: (params, ticket, callback) => {
          return _formatSampleTooltip(params.data[1]);
        },
        position: "right",
        showDelay: 60,
        hideDelay: 20,
      },
    },
  ];
  var _counts_values =
    DATA.meta_components.variant_counts_clustering.counts.map((_) => _[2]);
  _counts_values.sort((a, b) => a - b);
  var _counts_lim =
    _counts_values.slice(0, Math.round(_counts_values.length * 0.95)).pop() + 1;
  visualMap.push({
    type: "continuous",
    seriesIndex: _manifold_series_idx,
    min: 1,
    max: _counts_lim,
    range: [1, _counts_lim],
    color: colorMapOrRd,
    orient: "horizontal",
    top: 25,
    left: "3%",
    itemWidth: 10,
    text: ["> " + (_counts_lim - 1), "No. variants: 1 "],
    ...textOption,
  });

  return {
    title: [
      {
        text: "Sample's Variant Count by Feature",
        ...titleOption,
        left: "3%",
      },
      {
        text:
          "t-SNE of Sample Distances, Trustworthiness: " +
          parseFloat(
            DATA.meta_components.sample_distance_manifold.trustworthiness
          ).toFixed(2),
        ...titleOption,
        left: "53%",
      },
    ],
    legend: legend,
    grid: [
      {
        top: 55,
        bottom: "10%",
        left: "3%",
        width: "30%",
        containLabel: true,
        show: true,
      },
      {
        top: 55,
        bottom: "10%",
        left: "33%",
        width: "15%",
        containLabel: true,
        show: false,
      },
      {
        top: 55,
        bottom: "10%",
        left: "53%",
        right: "3%",
        containLabel: true,
        show: false,
      },
    ],
    toolbox: {
      top: "top",
      right: "right",
      feature: {
        saveAsImage: {
          title: "Save as image",
          type: "png",
          pixelRatio: 2,
        },
        restore: {
          title: "Restore",
        },
      },
    },
    xAxis: [
      {
        type: "category",
        name: "Feature",
        gridIndex: 0,
        data: _feature_names.sort(
          (a, b) =>
            DATA.meta_components.variant_counts_clustering.feature_labels[a] -
            DATA.meta_components.variant_counts_clustering.feature_labels[b]
        ),
        ...axisOption,
        splitLine: {
          show: false,
          interval: 0,
          lineStyle: {
            color: "#ECECEF",
          },
        },
      },
      {
        type: "value",
        name: "Count",
        gridIndex: 1,
        ...axisOption,
      },
      {
        type: "value",
        name: "tSNE 1",
        gridIndex: 2,
        ...axisOption,
      },
    ],
    yAxis: [
      {
        type: "category",
        name: "Sample",
        gridIndex: 0,
        data: _sample_names_filtered.sort(
          (a, b) =>
            DATA.meta_components.variant_counts_clustering.sample_labels[a] -
            DATA.meta_components.variant_counts_clustering.sample_labels[b]
        ),
        ...axisOption,
        splitLine: {
          show: false,
          interval: 0,
          lineStyle: {
            color: "#ECECEF",
          },
        },
        nameGap: 100,
      },
      {
        type: "category",
        gridIndex: 1,
        data: _sample_names_filtered.sort(
          (a, b) =>
            DATA.meta_components.variant_counts_clustering.sample_labels[a] -
            DATA.meta_components.variant_counts_clustering.sample_labels[b]
        ),
        axisLabel: {
          show: false,
        },
        show: false,
      },
      {
        type: "value",
        name: "tSNE 2",
        gridIndex: 2,
        ...axisOption,
      },
    ],
    visualMap: visualMap,
    tooltip: {
      ...tooltipOption,
    },
    dataZoom: [
      {
        type: "inside",
        yAxisIndex: [0, 1],
      },
      {
        type: "slider",
        xAxisIndex: [0],
        ...sliderOption,
        bottom: 10,
        height: 10,
      },
      {
        type: "inside",
        xAxisIndex: [2],
        filterMode: "empty",
      },
      {
        type: "inside",
        yAxisIndex: [2],
        filterMode: "empty",
      },
    ],
    series: [
      ...Object.values(_manifold_series),
      ...Object.values(_counts_series),
    ],
  };
}

/**
 * Constructs the ECharts option specification used to fill the dashboard for the features record type.
 *
 * @returns {Object} ECharts option specification to display via the `Dashboard` class defined in the `musial.web.results.js` module.
 */
function getDashboardOptionFeatures() {
  var _feature_names = Object.keys(DATA.storage.features);
  var _sample_names_filtered =
    SAMPLES_ACTIVE.length > 0
      ? SAMPLES_ACTIVE
      : Object.keys(DATA.storage.samples);
  var visualMap = []; // Stores the visual map specifications for the ECharts option.
  var legend = []; // Stores the legend specifications for the ECharts option.

  // Construct series for variants counts.
  var _counts_series = [
    {
      id: "dashboard.Features.counts",
      type: "heatmap",
      xAxisIndex: 0,
      yAxisIndex: 0,
      data: DATA.meta_components.variant_counts_clustering.counts
        .map((_) => [_[0], _[1], _[2]])
        .filter((_) => _[2] > 0),
      z: 0,
      tooltip: {
        formatter: (params, ticket, callback) => {
          let _ = _formatFeatureTooltip(params.data[1]);
          _ += _formatSampleTooltip(params.data[0], params.data[1]);
          return _;
        },
        position: "right",
        showDelay: 60,
        hideDelay: 20,
      },
      progressive: 2000,
      progressiveThreshold: 1000,
    },
    {
      id: "dashboard.Features.marginal",
      type: "bar",
      xAxisIndex: 1,
      yAxisIndex: 1,
      data: Object.entries(DATA.storage.features).map((_) => {
        let variant_count = 0;
        Object.values(_[1].nucleotideVariants).forEach((_) => {
          variant_count += Object.keys(_).length;
        });
        return [variant_count, _[0]];
      }),
      itemStyle: {
        color: "#333333",
      },
      tooltip: {
        formatter: (params, ticket, callback) => {
          return _formatFeatureTooltip(params.data[1]);
        },
        position: "right",
        showDelay: 60,
        hideDelay: 20,
      },
    },
  ];
  var _counts_values =
    DATA.meta_components.variant_counts_clustering.counts.map((_) => _[2]);
  _counts_values.sort((a, b) => a - b);
  var _counts_lim =
    _counts_values.slice(0, Math.round(_counts_values.length * 0.95)).pop() + 1;
  visualMap.push({
    type: "continuous",
    seriesIndex: 0,
    min: 1,
    max: _counts_lim,
    range: [1, _counts_lim],
    color: colorMapOrRd,
    orient: "horizontal",
    top: 25,
    left: "3%",
    itemWidth: 10,
    text: ["> " + (_counts_lim - 1), "No. variants: 1 "],
    ...textOption,
  });

  // Construct series for allele network.
  let _feature_name = $("#results-chart§control-features§feature")[0].value;
  let _allele_network_color_by = $("#results-chart§control-features§color")[0]
    .value;
  let _allele_network_categories = [];
  let _allele_network_categories_mapping = {};
  if (_allele_network_color_by == "table_filter") {
    _allele_network_categories.push({
      name: "Pass",
      itemStyle: { color: "#de3c4b" },
    });
    _allele_network_categories.push({
      name: "Other",
      itemStyle: { color: "#333333" },
    });
    _allele_network_categories_mapping["Pass"] = 0;
    _allele_network_categories_mapping["Other"] = 1;
  } else if (_allele_network_color_by == "proteoform") {
    if (DATA.storage.features[_feature_name].type == "non_coding") {
      _allele_network_categories.push({
        name: "N/A",
        itemStyle: { color: "#333333" },
      });
      _allele_network_categories_mapping["N/A"] = 0;
    } else {
      var _proteoform_colors = {
        get: (value) => {
          if (!_proteoform_colors._map.hasOwnProperty(value)) {
            _proteoform_colors._map[value] =
              _proteoform_colors._colors[_proteoform_colors._c];
            _proteoform_colors._c += 1;
            if (_proteoform_colors._c >= _proteoform_colors._colors.length)
              _proteoform_colors._c = 0;
          }
          return _proteoform_colors._map[value];
        },
        _map: {},
        _colors: [
          "#333333",
          "#a6cee3",
          "#1f78b4",
          "#b2df8a",
          "#33a02c",
          "#fb9a99",
          "#e31a1c",
          "#fdbf6f",
          "#ff7f00",
          "#cab2d6",
          "#6a3d9a",
        ],
        _c: 0,
      };
      let _idx = 0;
      Object.keys(DATA.storage.features[_feature_name].proteoforms).forEach(
        (proteoformName) => {
          let _clr = _proteoform_colors.get(proteoformName);
          _allele_network_categories.push({
            name: proteoformName,
            itemStyle: { color: _clr },
          });
          _allele_network_categories_mapping[proteoformName] = _idx;
          _idx++;
        }
      );
    }
  }
  legend.push({
    show: true,
    type: "scroll",
    top: 22,
    left: "53%",
    width: "40%",
    data: Object.keys(_allele_network_categories_mapping).map((_) => {
      return {
        name: _,
        icon: "circle",
      };
    }),
    itemGap: 25,
    ...textOption,
  });
  var _allele_network_series = {
    type: "graph",
    coordinateSystem: "cartesian2d",
    xAxisIndex: 2,
    yAxisIndex: 2,
    layout: "none",
    animation: false,
    roam: true,
    nodes: DATA.meta_components.allele_networks[_feature_name].nodes.map(
      (node) => {
        let _node_category;
        if (_allele_network_color_by == "table_filter") {
          _node_category = DATA.storage.features[_feature_name].alleles[
            node.name
          ].occurrence.some((_) => SAMPLES_ACTIVE.includes(_))
            ? _allele_network_categories_mapping["Pass"]
            : _allele_network_categories_mapping["Other"];
        } else {
          _node_category = _allele_network_categories_mapping[node._info[1]];
        }
        return {
          name: node.name,
          value: node.value,
          _info: node._info,
          category: _node_category,
          label:
            node.name == "reference"
              ? {
                  show: true,
                  position: "top",
                  distance: 1,
                  ...textOption.textStyle,
                }
              : {},
          symbolSize: node.name.startsWith("i")
            ? 0
            : Math.max(4, 3 * Math.log(node._info[0])),
        };
      }
    ),
    categories: _allele_network_categories,
    edges: DATA.meta_components.allele_networks[_feature_name].edges.map(
      (edge) => {
        return {
          source: edge.source,
          target: edge.target,
          value: edge.value,
          _info: edge._info,
          lineStyle: {
            color: "#AAAAAA",
          },
        };
      }
    ),
    edgeSymbol: ["circle", "none"],
    edgeSymbolSize: [5, 0],
    lineStyle: {
      curveness: 0.08,
    },
    emphasis: {
      disabled: true,
    },
    tooltip: {
      formatter: (params, ticket, callback) => {
        if (params.dataType == "node")
          return _formatAlleleTooltip(params.name, _feature_name);
        else if (params.dataType == "edge")
          return (
            `Allele <code>` +
            params.data._info[0] +
            `</code> and <code>` +
            params.data._info[1] +
            `</code> separated by <code>` +
            params.data.value +
            `</code> variant(s).</h6>`
          );
      },
      position: "left",
      showDelay: 60,
      hideDelay: 20,
    },
    legend: legend,
  };

  return {
    title: [
      {
        text: "Feature's Variant Count by Sample",
        ...titleOption,
        left: "3%",
      },
      {
        text: _feature_name + " Allele Network",
        ...titleOption,
        left: "53%",
      },
    ],
    legend: legend,
    grid: [
      {
        top: 55,
        bottom: "10%",
        left: "3%",
        width: "30%",
        containLabel: true,
        show: true,
      },
      {
        top: 55,
        bottom: "10%",
        left: "33%",
        width: "15%",
        containLabel: true,
        show: true,
      },
      {
        top: 55,
        bottom: "10%",
        left: "53%",
        right: "3%",
        containLabel: false,
        show: false,
      },
    ],
    toolbox: {
      top: "top",
      right: "right",
      feature: {
        saveAsImage: {
          title: "Save as image",
          type: "png",
          pixelRatio: 2,
        },
        restore: {
          title: "Restore",
        },
      },
    },
    tooltip: {},
    xAxis: [
      {
        type: "category",
        name: "Sample",
        gridIndex: 0,
        data: _sample_names_filtered.sort(
          (a, b) =>
            DATA.meta_components.variant_counts_clustering.sample_labels[a] -
            DATA.meta_components.variant_counts_clustering.sample_labels[b]
        ),
        inverse: true,
        ...axisOption,
        splitLine: {
          show: false,
          interval: 0,
          lineStyle: {
            color: "#ECECEF",
          },
        },
      },
      {
        type: "value",
        name: "Count",
        gridIndex: 1,
        ...axisOption,
      },
      {
        type: "value",
        gridIndex: 2,
        show: false,
      },
    ],
    yAxis: [
      {
        type: "category",
        name: "Feature",
        gridIndex: 0,
        data: _feature_names.sort(
          (a, b) =>
            DATA.meta_components.variant_counts_clustering.feature_labels[a] -
            DATA.meta_components.variant_counts_clustering.feature_labels[b]
        ),
        inverse: true,
        ...axisOption,
        splitLine: {
          show: false,
          interval: 0,
          lineStyle: {
            color: "#ECECEF",
          },
        },
        nameGap: 50,
      },
      {
        type: "category",
        gridIndex: 1,
        data: _feature_names.sort(
          (a, b) =>
            DATA.meta_components.variant_counts_clustering.feature_labels[a] -
            DATA.meta_components.variant_counts_clustering.feature_labels[b]
        ),
        axisLabel: {
          show: false,
        },
        inverse: true,
        show: false,
      },
      {
        type: "value",
        gridIndex: 2,
        show: false,
      },
    ],
    tooltip: {
      ...tooltipOption,
    },
    visualMap: visualMap,
    series: [..._counts_series, _allele_network_series],
    dataZoom: [
      {
        type: "inside",
        yAxisIndex: [0, 1],
      },
      {
        type: "slider",
        xAxisIndex: [0],
        ...sliderOption,
        bottom: 10,
        height: 10,
      },
      {
        type: "inside",
        xAxisIndex: [2],
        filterMode: "empty",
      },
      {
        type: "inside",
        yAxisIndex: [2],
        filterMode: "empty",
      },
    ],
  };
}

/**
 * Constructs the ECharts option specification used to fill the dashboard for the nucleotide variants record type.
 *
 * @returns {Object} ECharts option specification to display via the `Dashboard` class defined in the `musial.web.results.js` module.
 */
function getDashboardOptionNVariants() {
  let _series = {}; // Stores the actual series specifications for the ECharts option.
  let _effects = []; // Stores the effects of the nucleotide variants. These are used as series names to display in the legend.
  // Process features.
  Object.keys(DATA.storage.features).forEach((name) => {
    let _start = DATA.storage.features[name].start - 1;
    let _end = DATA.storage.features[name].end - 1;
    let _length = _end - _start + 1;
    _series[name] = {
      id: "dashboard.NVariants.featureTrack." + name,
      name: name,
      type: "line",
      xAxisIndex: 0,
      yAxisIndex: 0,
      data: [
        {
          value: [_start, 0],
          symbolSize: 0,
        },
        {
          value: [_start + Math.round((_end - _start) / 2), 0],
          symbolSize: 0,
          _zoomStart: _start - Math.round(_length * 0.02),
          _zoomEnd: _end + Math.round(_length * 0.02),
          label: {
            show: true,
            position: "inside",
            formatter: "{a}",
            align: "center",
            verticalAlign: "middle",
            color: "#333333",
            fontWeight: "bold",
          },
        },
        {
          value: [_end, 0],
          symbolSize: 0,
        },
      ],
      showSymbol: true,
      lineStyle: {
        color: "#cbd0e0",
        width: 20,
        type: "solid",
        shadowColor: "#111111",
        shadowBlur: 4,
        shadowOffsetX: 0,
        shadowOffsetY: 0,
        cap: "round",
        join: "miter",
      },
      tooltip: {
        formatter: (params, ticket, callback) => _formatFeatureTooltip(name),
      },
    };
  });
  // Process nucleotide variant records. Filter records by active samples.
  DATA.table_records.nucleotide_variants.rows
    .filter(
      (row) =>
        SAMPLES_ACTIVE.length == 0 ||
        row.occurrence.split(",").some((_) => SAMPLES_ACTIVE.includes(_))
    )
    .forEach((row) => {
      let alt = row.alternate_content;
      if (alt.length > 1) alt = "InDel";
      let effect = row.snpeff_Effect;
      if (effect == null) effect = "Unknown";
      else {
        effect = effect
          .split(new RegExp("_|&"))
          .map((_) => _.charAt(0).toUpperCase() + _.slice(1))
          .join(" ");
      }
      if (!_series.hasOwnProperty(effect)) {
        _series[effect] = {
          id: "dashboardNVariants." + effect,
          name: effect,
          type: "bar",
          xAxisIndex: 1,
          yAxisIndex: 1,
          barMinWidth: 2,
          barWidth: 5,
          barMaxWidth: 20,
          barGap: "1%",
          data: [],
          showBackground: true,
          backgroundStyle: {
            color: "#eff0f8",
          },
          stack: "position",
        };
        _effects.push(effect);
      }
      _series[effect].data.push({
        value: [row.position, row.frequency],
        itemStyle: {
          color: _mapNucleotideColor(alt),
        },
        tooltip: {
          formatter: (params, ticket, callback) => _formatNVariantTooltip(row),
        },
      });
    });
  // Construct the ECharts option specification; Incoorporate the series data constructed above.
  return {
    title: [
      {
        text: "Nucleotide Variants by Genome Position",
        ...titleOption,
      },
    ],
    grid: [
      {
        top: 50,
        height: 50,
        left: 50,
        right: 50,
      },
      {
        top: 105,
        bottom: 50,
        left: 50,
        right: 50,
      },
    ],
    legend: [
      {
        show: true,
        type: "scroll",
        top: 25,
        left: "center",
        width: "75%",
        data: _effects.map((effect) => {
          return {
            name: effect,
          };
        }),
        icon: "none",
        itemStyle: { color: "#333333" },
        itemGap: 25,
        itemWidth: 0,
        selector: [
          {
            type: "all",
            title: "Select all",
          },
          {
            type: "inverse",
            title: "Invert selection",
          },
        ],
        selectorLabel: {
          ...textOption.textStyle,
        },
        ...textOption,
      },
    ],
    toolbox: {
      top: "top",
      right: "right",
      feature: {
        saveAsImage: {
          title: "Save as image",
          type: "png",
          pixelRatio: 2,
        },
        restore: {
          title: "Restore",
        },
      },
    },
    xAxis: [
      {
        type: "value",
        gridIndex: 0,
        min: 1,
        max: parseInt(DATA.storage.info.reference_length),
        minInterval: 1,
        show: false,
      },
      {
        type: "value",
        gridIndex: 1,
        min: 1,
        max: parseInt(DATA.storage.info.reference_length),
        minInterval: 1,
        name: "Genome Position",
        ...axisOption,
        axisPointer: {
          show: true,
          snap: true,
          triggerEmphasis: false,
          triggerTooltip: false,
          ...labelOption,
        },
      },
    ],
    yAxis: [
      {
        type: "category",
        gridIndex: 0,
        data: [0],
        name: "Feature",
        ...axisOption,
        axisLabel: {
          show: false,
        },
        axisTick: {
          show: false,
        },
      },
      {
        type: "value",
        gridIndex: 1,
        name: "Frequency in Samples (%)",
        ...axisOption,
      },
    ],
    tooltip: {
      ...tooltipOption,
    },
    series: Object.values(_series),
    dataZoom: [
      {
        type: "inside",
        xAxisIndex: 0,
        filterMode: "none",
      },
      {
        type: "inside",
        xAxisIndex: [0, 1],
        filterMode: "weakFilter",
      },
      {
        type: "slider",
        yAxisIndex: [1],
        filterMode: "none",
        ...sliderOption,
        width: 10,
        right: 10,
      },
    ],
  };
}

/**
 * Constructs the ECharts option specification used to fill the dashboard for the amino acid variants record type.
 *
 * @returns {Object} ECharts option specification to display via the `Dashboard` class defined in the `musial.web.results.js` module.
 */
function getDashboardOptionAVariants() {
  return { };
}

/**
 * Maps nucleotide bases to HEX color codes
 *
 * @param {string} n Nuclotide base. One of "A", "C", "G", "T", "N", or "InDel".
 * @returns A HEX color code for the given nucleotide base.
 */
function _mapNucleotideColor(n) {
  if (n == "A") return "#1f78b4";
  else if (n == "C") return "#b2df8a";
  else if (n == "G") return "#a6cee3";
  else if (n == "T") return "#33a02c";
  else if (n == "InDel") return "#984ea3";
  else return "#333333";
}

function _formatSampleTooltip(name, feature) {
  let body = `<h6><code>Sample ` + name + `</code></h6>`;
  body += `<table>`;
  if (feature != undefined) {
    body +=
      `<tr><td style="padding-right: 6px">Allele</td><td style="padding-right: 6px">` +
      DATA.storage.samples[name].allele[feature] +
      `</td></tr>`;
    body +=
      `<tr><td style="padding-right: 6px">Proteoform</td><td style="padding-right: 6px">` +
      DATA.storage.samples[name].proteoform[feature] +
      `</td></tr>`;
  }
  for (const [k, v] of Object.entries(DATA.storage.samples[name].info)) {
    if (v == null) continue;
    let _k = k;
    _k = _k.replace(/_/g, " ");
    _k = _k.charAt(0).toUpperCase() + _k.slice(1);
    body +=
      `<tr><td style="padding-right: 6px">` +
      _k +
      `</td><td style="padding-right: 6px">` +
      v +
      `</td></tr>`;
  }
  body += `</table>`;
  return body;
}

/**
 * Formats the tooltip for a feature.
 *
 * @param {string} name The name of the feature.
 * @returns HTML string; formatted tooltip.
 */
function _formatFeatureTooltip(name) {
  let body = `<h6><code>Feature ` + name + `</code></h6>`;
  body += `<table>`;
  body +=
    `<tr><td style="padding-right: 6px">Start</td><td style="padding-right: 6px">` +
    DATA.storage.features[name].start +
    `</td></tr>`;
  body +=
    `<tr><td style="padding-right: 6px">End</td><td style="padding-right: 6px">` +
    DATA.storage.features[name].end +
    `</td></tr>`;
  body +=
    `<tr><td style="padding-right: 6px">Strand</td><td style="padding-right: 6px">` +
    (DATA.storage.features[name].isSense ? "+" : "-") +
    `</td></tr>`;
  body +=
    `<tr><td style="padding-right: 6px">No. Alleles</td><td style="padding-right: 6px">` +
    Object.keys(DATA.storage.features[name].alleles).length +
    `</td></tr>`;
  if (DATA.storage.features[name].hasOwnProperty("proteoforms"))
    body +=
      `<tr><td style="padding-right: 6px">No. Proteoforms</td><td style="padding-right: 6px">` +
      Object.keys(DATA.storage.features[name].proteoforms).length +
      `</td></tr>`;
  for (const [k, v] of Object.entries(DATA.storage.features[name].info)) {
    if (v == null) continue;
    if (k.startsWith("number_of_")) continue;
    let _k = k;
    _k = _k.replace(/_/g, " ");
    _k = _k.charAt(0).toUpperCase() + _k.slice(1);
    body +=
      `<tr><td style="padding-right: 6px">` +
      _k +
      `</td><td style="padding-right: 6px">` +
      v +
      `</td></tr>`;
  }
  body += `</table>`;
  return body;
}

function _formatAlleleTooltip(name, feature) {
  let body = `<h6><code>Allele ` + name + `</code></h6>`;
  body += `<table>`;
  for (const [k, v] of Object.entries(
    DATA.storage.features[feature].alleles[name].info
  )) {
    if (v == null) continue;
    if (k.startsWith("number_of_")) continue;
    let _k = k;
    _k = _k.replace(/_/g, " ");
    _k = _k.charAt(0).toUpperCase() + _k.slice(1);
    body +=
      `<tr><td style="padding-right: 6px">` +
      _k +
      `</td><td style="padding-right: 6px">` +
      v +
      `</td></tr>`;
  }
  let allele_variants =
    DATA.storage.features[feature].alleles[name].variants.split(",");
  if (allele_variants.length > 10)
    allele_variants = allele_variants.slice(0, 10).concat(["..."]);
  body +=
    `<tr><td style="padding-right: 6px">Variants</td><td style="padding-right: 6px">` +
    allele_variants.join(" ") +
    `</td></tr>`;
  body += `</table>`;
  return body;
}

/**
 * Formats the tooltip for a nucleotide variant record.
 *
 * @param {Object} record The nucleotide variant record.
 * @returns HTML string; formatted tooltip.
 */
function _formatNVariantTooltip(record) {
  let body =
    `<h6><code>` +
    record.alternate_content +
    ` at ` +
    record.position +
    `</code></h6>`;
  body += `<table>`;
  for (const [k, v] of Object.entries(record)) {
    if (v == null) continue;
    if (
      k.startsWith("alternate_content") ||
      k.startsWith("occurrence") ||
      k.startsWith("position")
    )
      continue;
    let _k = k;
    _k = _k.replace(/_/g, " ");
    _k = _k.charAt(0).toUpperCase() + _k.slice(1);
    body +=
      `<tr><td style="padding-right: 6px">` +
      _k +
      `</td><td style="padding-right: 6px">` +
      v +
      `</td></tr>`;
  }
  body += `</table>`;
  return body;
}

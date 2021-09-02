let viewHeight = 0;
var jsonInput;
var echart;
var echartOption;
var labelsOn = false;
var viewer;
var viewerTooltipLabel;
var showLabels = {
    valueInternal: false,
    valueListener: function(val) {},
    set value(val) {
        this.valueInternal = val;
        this.valueListener(val);
    },
    get value() {
        return this.valueInternal;
    },
    registerListener: function(listener) {
        this.valueListener = listener;
    }
};
showLabels.registerListener( function( val ) {
    if ( val ) {
        echart.setOption( {
            series: [
                {
                    name: 'Reference',
                    label: {
                        show: true
                    }
                },{
                    name: 'Protein',
                    label: {
                        show: true
                    }
                }
            ]
        } );
    } else {
        echart.setOption( {
            series: [
                {
                    name: 'Reference',
                    label: {
                        show: false
                    }
                },{
                    name: 'Protein',
                    label: {
                        show: false
                    }
                }
            ]
        } );
    }
});

function fileInputHandler( event ) {
    document.getElementById( "fileInput" ).style.display = "none";
    document.getElementById( "fileInputIcon" ).style.display = "none";
    document.getElementById( "fileInputText" ).style.display = "none";
    document.getElementById( "fileInputProcessLoader" ).style.display = "flex";

    var reader = new FileReader();
    reader.onload = function(event) {
        jsonInput = JSON.parse(event.target.result);
        processJSONInput( jsonInput );
        document.getElementById( "input" ).style.display = "none";
        document.getElementById( "view" ).style.display = "flex";
    }
    reader.readAsText(event.target.files[0]);
}

function processJSONInput( jsonObj ) {
    if ( 'ProteinStructure' in jsonObj ) {
        loadProteinView( );
    } else {
        document.getElementById( "proteinView" ).style.display = "none";
    }
    loadEChartsView( );
}

window.onload = _ => {
    new ResizeObserver( ( entries ) => {
        viewHeight = entries[ 0 ].contentRect.height;
    } ).observe( document.getElementById( "view" ) );

    new ResizeObserver( ( entries ) => {
        document.getElementById( "eChartsView" ).style.height = ( viewHeight - entries[ 0 ].contentRect.height ).toString( ) + "px";
        echart.resize( );
    } ).observe( document.getElementById( "proteinView" ) );
};

function dnaIntToChar( i ) {
    let dnaSymbols = "ACGT$.:;-aqycvdgfhtrz".split( "" );
    if ( i > 0 & i < 21 ) {
        let dnaSymbol = dnaSymbols[ i - 1 ];
        if ( "Aaqy".split( "" ).includes( dnaSymbol ) ) {
            return "A";
        } else if ( "Ccvd".split( "" ).includes( dnaSymbol ) ) {
            return "C";
        } else if ( "Ggfh".split( "" ).includes( dnaSymbol ) ) {
            return "G";
        } else if ( "Ttrz".split( "" ).includes( dnaSymbol ) ) {
            return "T";
        } else if ( ".:;".split( "" ).includes( dnaSymbol ) ) {
            return "Reference";
        } else if ( dnaSymbol == "-" ) {
            return "Deletion";
        } else if ( dnaSymbol == "$" ) {
            return "No information";
        }
    } else {
        return "";
    }
}

function dnaCharToInt( c ) {
    let dnaSymbols = "ACGT$.:;-aqycvdgfhtrz".split( "" );
    return ( dnaSymbols.indexOf( c ) + 1 );
}

function aaIntToChar( i ) {
    let aaSymbols = "HKRDESTNQAVLIMFWYPGC$&-".split( "" );
    if ( i > 0 & i < 24 ) {
        let aaSymbol = aaSymbols[ i - 1 ];
        if ( aaSymbol == "-" ) {
            return "Deletion";
        } else if ( "$&".split( "" ).includes( aaSymbol ) ) {
            return "No information";
        } else {
            return aaSymbols[ i - 1 ];
        }
    } else {
        return "";
    }
}

function aaCharToInt( c ) {
    let aaSymbols = "HKRDESTNQAVLIMFWYPGC$&-".split( "" );
    return ( aaSymbols.indexOf( c ) + 1 );
}

function loadEChartsView( ) { 
    var chartDom = document.getElementById('eChartsView');
    echart = echarts.init(chartDom, { "renderer": "svg" } );

    var referenceLabel = [ "Ref. Gene" ];
    let referenceSequence = jsonInput.ReferenceSequence;
    var length = referenceSequence.length;
    var referenceData = Array.from(Array(length).keys()).map( position => [ position, 0, dnaCharToInt( referenceSequence[ position ] ) ] );

    var positions = Array.from(Array(length).keys()).map( position => ( position + 1 ).toString( ) );

    var proteinLabel = [ "Ref. Protein" ];
    let proteinSequence = jsonInput.ProteinSequence;
    let frameShiftStart = 0;
    let frameShiftEnd = 0;
    if ( proteinSequence.startsWith( "&&" ) ) {
        proteinSequence = proteinSequence.substring( 2, proteinSequence.length );
        frameShiftStart = 2;
    } else if ( proteinSequence.startsWith( "&" ) ) {
        proteinSequence = proteinSequence.substring( 1, proteinSequence.length );
        frameShiftStart = 1;
    }
    if ( proteinSequence.endsWith( "&&" ) ) {
        proteinSequence = proteinSequence.substring( 0, proteinSequence.length - 2 );
        frameShiftEnd = 2;
    } else if ( proteinSequence.endsWith( "&" ) ) {
        proteinSequence = proteinSequence.substring( 0, proteinSequence.length - 1 );
        frameShiftEnd = 1;
    }
    var proteinSequenceData = Array.from(Array(length - frameShiftStart - frameShiftEnd).keys()).map( position => [ position + frameShiftStart, 0, aaCharToInt( proteinSequence[ Math.ceil( ( position + 1 ) / 3 ) - 1 ] ) ] );
    for ( let i = 0; i < frameShiftStart; i ++ ) {
        proteinSequenceData.unshift( [ i, 0, aaCharToInt( "&" ) ] );
    }
    for ( let i = 0; i < frameShiftEnd; i ++ ) {
        proteinSequenceData.push( [ length + i, 0, aaCharToInt( "&" ) ] );
    }

    var samplesLabels = jsonInput.SampleNames;
    var samplesData = [ ];
    var samplesAnnotation = [ ];
    var snvCounts = {
        "Ref.": [ ],
        "Ref. discarded": [ ],
        "Del.": [ ],
        "A": [ ],
        "A discarded": [ ],
        "C": [ ],
        "C discarded": [ ],
        "G": [ ],
        "G discarded": [ ],
        "T": [ ],
        "T discarded": [ ]
    };
    var variantPositions = Object.keys( jsonInput.PerPositionVariants );
    for ( let j = 0; j < length; j ++ ) {
        let variantPosition = false;
        if ( variantPositions.includes( String( j + 1 ) ) ) {
            variantPosition = true;
        }
        for ( let i = 0; i < samplesLabels.length; i++ ) {
            let dnaChar = ".";
            let annotation = "";
            if ( variantPosition ) {
                dnaChar = jsonInput.PerPositionVariants[ String( j + 1 ) ][ i ];
                annotation = jsonInput.PerPositionAnnotations[ String( j + 1 ) ][ i ];
                samplesData.push( [ j, i, dnaCharToInt( dnaChar ) ] );
                samplesAnnotation.push( [ j, i, annotation ] );
            }    
        }
        if ( variantPosition ) {
            snvCounts[ "Ref." ].push( parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 0 ] ) );
            snvCounts[ "Ref. discarded" ].push( 
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 1 ] ) + 
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 2 ] )
            );
            snvCounts[ "Del." ].push( parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 3 ] ) );
            snvCounts[ "A" ].push( parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 4 ] ) );
            snvCounts[ "A discarded" ].push(
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 5 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 6 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 7 ] )
            );
            snvCounts[ "C" ].push( parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 8 ] ) );
            snvCounts[ "C discarded" ].push(
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 9 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 10 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 11 ] )
            );
            snvCounts[ "G" ].push( parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 12 ] ) );
            snvCounts[ "G discarded" ].push(
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 13 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 14 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 15 ] )
            );
            snvCounts[ "T" ].push( parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 16 ] ) );
            snvCounts[ "T discarded" ].push(
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 17 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 18 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( j + 1 ) ][ 19 ] )
            );
        } else {
            Object.keys( snvCounts ).forEach( ( key ) => {
                snvCounts[ key ].push( null );
            } );
        }
    }

    echartOption = {
        animation: false,
        tooltip: {
            confine: true,
            trigger: 'item',
            triggerOn: 'mousemove',
            alwaysShowContent: true,
            position: [ 10, 10 ],
            extraCssText: 'width: 10%; height: 92%',
            textStyle: {
                fontSize: 11
            },
            transitionDuration: 3,
            borderWidth: 2,
            borderColor: "#E4E5ED",
            formatter: function (params, ticket, callback) {
                let seriesIndex = params.seriesIndex;
                let tooltipHtmlString = "";
                let position = 0;
                let referenceContent = "";
                let proteinSequenceContent = "";
                let sampleName = "-";
                let sampleContent = "-";
                let sampleAnnotation = [];
                let aCount = 0;
                let cCount = 0;
                let gCount = 0;
                let tCount = 0;
                let delCount = 0;
                let refCount = 0;
                let discardedCount = 0;
                let dataIndex = 0;
                if ( [2,3,4,5,6,7,8,9,10,11,12].includes( seriesIndex ) ) {
                    // CASE: Item is from counts chart.
                    dataIndex = params.dataIndex;
                } else {
                    // CASE: Item is from other chart.
                    dataIndex = params.data[ 0 ];
                    
                }
                position = dataIndex + 1;
                referenceContent = dnaIntToChar( referenceData[ dataIndex ][ 2 ] );
                proteinSequenceContent = aaIntToChar( proteinSequenceData[ dataIndex ][ 2 ] );
                if ( seriesIndex == 1 ) {
                    // CASE: Item is from samples chart.
                    let index = params.data[ 1 ];
                    sampleName = samplesLabels[ index ];
                    sampleContent = dnaIntToChar( params.data[ 2 ] );
                    let annotationString = samplesAnnotation[ dataIndex ][ 2 ];
                    if ( annotationString != "." ) {
                        for ( let annotation of annotationString.split( ";" ) ) {
                            let annotationFields = annotation.split( "=" );
                            sampleAnnotation.push( annotationFields );
                        }
                    }
                }
                if ( variantPositions.includes( String( position ) ) ) {
                    aCount = snvCounts[ "A" ][ dataIndex ];
                    cCount = snvCounts[ "C" ][ dataIndex ];
                    gCount = snvCounts[ "G" ][ dataIndex ];
                    tCount = snvCounts[ "T" ][ dataIndex ];
                    delCount = snvCounts[ "Del." ][ dataIndex ];
                    refCount = snvCounts[ "Ref." ][ dataIndex ];
                    discardedCount = (
                        snvCounts[ "A discarded" ][ dataIndex ] +
                        snvCounts[ "C discarded" ][ dataIndex ] +
                        snvCounts[ "G discarded" ][ dataIndex ] +
                        snvCounts[ "T discarded" ][ dataIndex ] +
                        snvCounts[ "Ref. discarded" ][ dataIndex ]
                    );  
                } else {
                    refCount = samplesLabels.length;
                }
                tooltipHtmlString += "<b>Reference Information</b><br>"
                tooltipHtmlString += "Position: " + position + "<br>";
                tooltipHtmlString += "Nucleotide: " + referenceContent + "<br>";
                tooltipHtmlString += "Amino Acid: " + proteinSequenceContent + "<br>";
                tooltipHtmlString += "<hr>";
                tooltipHtmlString += "<b>Sample Information</b><br>"
                tooltipHtmlString += "Sample Name: " + sampleName + "<br>";
                tooltipHtmlString += "Sample Variant: " + sampleContent + "<br>";
                for ( let sa of sampleAnnotation ) {
                    tooltipHtmlString += sa[ 0 ] + ": " + sa[ 1 ] + "<br>";
                }
                tooltipHtmlString += "<hr>";
                tooltipHtmlString += "<b>Counts</b><br>"
                tooltipHtmlString += "# A: " + aCount + "<br>";
                tooltipHtmlString += "# C: " + cCount + "<br>";
                tooltipHtmlString += "# G: " + gCount + "<br>";
                tooltipHtmlString += "# T: " + tCount + "<br>";
                tooltipHtmlString += "# Ref.: " + refCount + "<br>";
                tooltipHtmlString += "# Del.: " + delCount + "<br>";
                tooltipHtmlString += "# Discarded: " + discardedCount;
                return tooltipHtmlString;
            }
        },
        /*
        toolbox: {
            orient: 'vertical',
            feature: {
                saveAsImage: {},
                dataView: {}
            }
        },
        */
        grid: [
            {
                top: '10',
                left: '17%',
                height: '5%',
                width: '77%'
            },
            {
                top: '50', //'14.5%',
                left: '17%',
                height: '5%',
                width: '77%'
            },
            {
                top: '90', //'21.5%',
                left: '17%',
                height: '54%',
                width: '77%',
                show: true,
                backgroundColor: '#EEEEEE'
            },
            {
                top: '410', // 73.5%',
                left: '17%',
                height: '16%',
                width: '77%',
                show: true,
                backgroundColor: '#EEEEEE'
            }
        ],
        xAxis: [
            {
                type: 'category',
                data: positions,
                splitArea: {
                    show: false
                },
                show: false,
                gridIndex: 1
            },
            {
                type: 'category',
                data: positions,
                splitArea: {
                    show: false
                },
                show: false,
                gridIndex: 2
            },
            {
                type: 'category',
                data: positions,
                splitArea: {
                    show: false
                },
                gridIndex: 3
            },
            {
                type: 'category',
                data: positions,
                splitArea: {
                    show: false
                },
                show: false,
                gridIndex: 0
            }
        ],
        yAxis: [
            {
                type: 'category',
                data: referenceLabel,
                splitArea: {
                    show: false
                },
                gridIndex: 1
            },
            {
                type: 'category',
                data: samplesLabels,
                axisTick: {
                    alignWithLabel: true
                },
                splitArea: {
                    show: false
                },
                gridIndex: 2,
                inverse: true
            },
            {
                type: 'value',
                min: 0,
                max: samplesLabels.length,
                splitArea: {
                    show: false
                },
                gridIndex: 3
            },
            {
                type: 'category',
                data: proteinLabel,
                splitArea: {
                    show: false
                },
                gridIndex: 0
            }
        ],
        visualMap: [
            {
                type: 'piecewise',
                pieces: [
                    // ACGT$.:;-aqygfhtrzcvd
                    { min: 1, max: 1, color: "#12e049" }, // A
                    { min: 2, max: 2, color: "#4B90DE" }, // C
                    { min: 3, max: 3, color: "#DEC54B" }, // G
                    { min: 4, max: 4, color: "#DF4B4B" }, // T
                    { min: 5, max: 5, color: "#ed4ce5" }, // No match with protein sequence.
                    { min: 6, max: 6, color: "#DDDDDD" }, // Reference
                    { min: 7, max: 7, color: "#999999" }, // Reference, low quality
                    { min: 8, max: 8, color: "#999999"}, // Reference, low coverage
                    { min: 9, max: 9, color: "#333333" }, // Deletion
                    { min: 10, max: 10, color: "#8bd6a0" }, // A, low quality
                    { min: 11, max: 11, color: "#8bd6a0" }, // A, low coverage
                    { min: 12, max: 12, color: "#8bd6a0" }, // A, low frequency
                    { min: 13, max: 13, color: "#a7bcd4" }, // C, low quality
                    { min: 14, max: 14, color: "#a7bcd4" }, // C, low coverage
                    { min: 15, max: 15, color: "#a7bcd4" }, // C, low frequency
                    { min: 16, max: 16, color: "#e0d59f" }, // G, low quality
                    { min: 17, max: 17, color: "#e0d59f" }, // G, low coverage
                    { min: 18, max: 18, color: "#e0d59f" }, // G, low frequency
                    { min: 19, max: 19, color: "#e39f9f" }, // T, low quality
                    { min: 20, max: 20, color: "#e39f9f" }, // T, low coverage
                    { min: 21, max: 21, color: "#e39f9f" } // T, low frequency
                ],
                seriesIndex: [ 0, 1 ],
                show: false
            },
            {
            type: 'piecewise',
                pieces: [
                    { min: 1, max: 3, color: "#ADF5FF" },
                    { min: 4, max: 5, color: "#D55672" },
                    { min: 6, max: 9, color: "#81D488" },
                    { min: 10, max: 14, color: "#EDC79B" },
                    { min: 15, max: 17, color: "#E5BEED" },
                    { min: 18, max: 19, color: "#93827F" },
                    { min: 20, max: 20, color: "#FFE66D" },
                    { min: 21, max: 22, color: "#ed4ce5" },
                    { min: 23, max: 23, color: "#333333" }
                ],
                seriesIndex: [ 13 ],
                show: false
            }
        ],
        legend: {
            show: false
        },
        dataZoom: [
            {
                type: 'slider',
                xAxisIndex: [ 0, 1, 3 ],
                minValueSpan: 30,
                maxValueSpan: length,
                realtime: false,
                throttle: 100,
                bottom: 10,
                id: "positionZoom",
                showDataShadow: false,
                backgroundColor: "transparent",
                fillerColor: "transparent",
                borderColor: "rgba(228, 229, 237, 0.733)",
                handleStyle: {
                    borderColor: "rgba(96, 113, 150, 1.0)"
                },
                moveHandleStyle: {
                    color: "rgba(96, 113, 150, 0.533)"
                },
                emphasis: {
                    moveHandleStyle: {
                        color: "rgba(96, 113, 150, 1.0)"
                    }
                }
            },
            {
                type: 'slider',
                yAxisIndex: [ 1 ],
                minValueSpan: 8,
                maxValueSpan: samplesLabels.length,
                throttle: 5,
                right: "4%",
                backgroundColor: "transparent",
                fillerColor: "transparent",
                borderColor: "rgba(228, 229, 237, 0.733)",
                handleStyle: {
                    borderColor: "rgba(96, 113, 150, 1.0)"
                },
                moveHandleStyle: {
                    color: "rgba(96, 113, 150, 0.533)"
                },
                emphasis: {
                    moveHandleStyle: {
                        color: "rgba(96, 113, 150, 1.0)"
                    }
                }
            }
        ],
        series: [
            {
                name: 'Reference',
                type: 'heatmap',
                data: referenceData,
                label: {
                    show: false,
                    formatter: function ( params ) {
                        return dnaIntToChar( params.value[ 2 ] );
                    },
                    color: "black"
                },
                emphasis: {
                    itemStyle: {
                        shadowBlur: 10,
                        shadowColor: 'rgba(0, 0, 0, 0.5)'
                    }
                },
                xAxisIndex: 0,
                yAxisIndex: 0,
                itemStyle: {
                    borderColor: "#E8E9ED",
                    borderWidth: 0.1
                }
            },
            {
                name: 'Samples',
                type: 'heatmap',
                data: samplesData,
                label: {
                    show: false,
                    formatter: function ( params ) {
                        if ( params.value[ 2 ] == 5 ) {
                            return "";
                        } else {
                            return dnaIntToChar( params.value[ 2 ] );
                        }
                    },
                    color: "black"
                },
                xAxisIndex: 1,
                yAxisIndex: 1,
                itemStyle: {
                    borderColor: "#E8E9ED",
                    borderWidth: 0.1
                },
                animation: false
            },
            {
                name: 'A calls',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#12e049"
                },
                data: snvCounts[ "A" ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'A discarded',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#8bd6a0"
                },
                data: snvCounts[ "A discarded" ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'C calls',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#4B90DE"
                },
                data: snvCounts[ "C" ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'C discarded',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#a7bcd4"
                },
                data: snvCounts[ "C discarded" ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'G calls',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#DEC54B"
                },
                data: snvCounts[ "G" ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'G discarded',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#e0d59f"
                },
                data: snvCounts[ "G discarded" ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'T calls',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#DF4B4B"
                },
                data: snvCounts[ "T" ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'T discarded',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#e39f9f"
                },
                data: snvCounts[ "T discarded" ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'Deletions',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#333333"
                },
                data: snvCounts[ "Del." ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'Reference discarded',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#999999"
                },
                data: snvCounts[ "Ref. discarded" ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'Reference calls',
                type: 'bar',
                stack: 'SNV_AGGREGATED',
                itemStyle: {
                    color: "#DDDDDD"
                },
                data: snvCounts[ "Ref." ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'Protein',
                type: 'heatmap',
                data: proteinSequenceData,
                label: {
                    show: false,
                    formatter: function ( params ) {
                        if ( params.value[ 0 ] % 3 == 1 ) {
                            return aaIntToChar( params.value[ 2 ] );
                        } else {
                            return "";
                        }
                    },
                    color: "black"
                },
                itemStyle: {
                    borderColor: "#E8E9ED",
                    borderWidth: 0.1
                },
                xAxisIndex: 3,
                yAxisIndex: 3
            }
        ]
    };

    echart.on('click', function (params) {
        let seriesIndex = params.seriesIndex;
        let dataIndex = undefined;
        if ( [ 0, 1, 13 ].includes( seriesIndex )  ) {
            dataIndex = params.data[ 0 ];
        } else {
            dataIndex = params.dataIndex;
        }
        let position = Math.ceil( ( dataIndex + 1 ) / 3 );
        viewer.removeAllLabels();
        viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
        viewer.addStyle({resi: [position]},{stick:{radius:1,color:"#FFBE0B"}});
        viewer.addResLabels({resi: [position]});
        viewer.render();
    });
    
    echart.on('datazoom', function(params) {
        if ( params.dataZoomId == "positionZoom" ) {
            startPosition = ( Math.round( params.start ) / 100 ) * length;
            endPosition = ( Math.round( params.end ) / 100 ) * length;
        }
        if ( ( endPosition - startPosition ) <= 100 ) {
            showLabels.value = true;
        } else {
            showLabels.value = false;
        }
    });

    echart.setOption(echartOption);
}

function loadProteinView( ) {
    let element = $('#proteinView');
    let config = { backgroundColor: '#FFFFFF' };
    viewer = $3Dmol.createViewer( element, config );
    viewer.addModel( jsonInput.ProteinStructure, "pdb" );
    viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
    viewer.setViewStyle({style:"outline"});
    viewer.zoomTo();
    viewer.render();
    viewer.zoom(1.5, 1000);
}
let viewHeight = 0;
var jsonInput;
var echart;
var echartOption;
var labelsOn = false;
var viewer;
var viewerTooltipLabel;
var proteinResidueNumberMap = { };
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
var colorSchemes = {
    "nucleotides": {
        "A": "#3E885B",
        "dA": "#93CDAA", 
        "C": "#0471A6",
        "dC": "#79D0FC",
        "G": "#E8C547",
        "dG": "#F3E1A0",
        "T": "#DB5461",
        "dT": "#ECA7AE",
        "del": "transparent",
        "ref": "#EEEEEE",
        "dref": "#9093A7",
        "nomatch": "transparent"
    },
    "aminoacids": {
        "H": "#87ABFF", // Polar (positive), Basic
        "K": "#87ABFF",
        "R": "#87ABFF",
        "D": "#F75050", // Polar (negative), Acidic
        "E": "#F75050",
        "S": "#AEDEAA", // Polar (neutral)
        "T": "#AEDEAA",
        "N": "#AEDEAA",
        "Q": "#AEDEAA",
        "C": "#AEDEAA",
        "F": "#FFEAB1", // Aromatic
        "W": "#FFEAB1",
        "Y": "#FFEAB1",
        "A": "#D387F3", // Aliphatic
        "V": "#D387F3",
        "L": "#D387F3",
        "I": "#D387F3",
        "M": "#D387F3",
        "P": "#D387F3",
        "G": "#D387F3",
        "del": "transparent",
        "nomatch": "transparent"
    },
    "mol3d": {
        "HIS": "#87ABFF", // Polar (positive), Basic
        "LYS": "#87ABFF",
        "ARG": "#87ABFF",
        "ASP": "#F75050", // Polar (negative), Acidic
        "GLU": "#F75050",
        "SER": "#AEDEAA", // Polar (neutral)
        "THR": "#AEDEAA",
        "ASN": "#AEDEAA",
        "GLN": "#AEDEAA",
        "CYS": "#AEDEAA",
        "PHE": "#FFEAB1", // Aromatic
        "TRP": "#FFEAB1",
        "TYR": "#FFEAB1",
        "ALA": "#D387F3", // Aliphatic
        "VAL": "#D387F3",
        "LEU": "#D387F3",
        "ILE": "#D387F3",
        "MET": "#D387F3",
        "PRO": "#D387F3",
        "GLY": "#D387F3"
    }
};

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
    let dnaSymbols = "ACGT?.:;-aqycvdgfhtrz".split( "" );
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
            return "-";
        } else if ( dnaSymbol == "?" ) {
            return "?";
        }
    } else {
        return "";
    }
}

function dnaCharToInt( c ) {
    let dnaSymbols = "ACGT?.:;-aqycvdgfhtrz".split( "" );
    return ( dnaSymbols.indexOf( c ) + 1 );
}

function aaIntToChar( i ) {
    let aaSymbols = "HKRDESTNQCFWYAVLIMPG-?".split( "" );
    if ( i > 0 & i < 23 ) {
        let aaSymbol = aaSymbols[ i - 1 ];
        if ( aaSymbol == "-" ) {
            return "-";
        } else if ( aaSymbol == "?" ) {
            return "?";
        } else {
            return aaSymbols[ i - 1 ];
        }
    } else {
        return "";
    }
}

function aaCharToInt( c ) {
    let aaSymbols = "HKRDESTNQCFWYAVLIMPG-?".split( "" );
    return ( aaSymbols.indexOf( c ) + 1 );
}

function loadEChartsView( ) { 
    var chartDom = document.getElementById('eChartsView');
    echart = echarts.init(chartDom, { "renderer": "svg" } );

    var referenceLabel = [ jsonInput.ReferenceName ];
    let referenceSequence = jsonInput.ReferenceSequence;
    var length = referenceSequence.length;
    var referenceData = Array.from(Array(length).keys()).map( position => [ position, 0, dnaCharToInt( referenceSequence[ position ] ) ] );

    var positions = Array.from(Array(length).keys()).map( position => ( position + 1 ).toString( ) );

    var proteinLabel = [ jsonInput.ProteinName ];
    let proteinSequence = jsonInput.ProteinSequence;
    var proteinData = Array.from(Array(length).keys()).map( position => [ position, 0, aaCharToInt( proteinSequence[ position ] ) ] );
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
        if ( variantPositions.includes( String( j + 1 ) ) ) variantPosition = true;
        if ( variantPosition ) {
            let index = j + 1;
            
            for ( let i = 0; i < samplesLabels.length; i++ ) {
                let dnaChar = jsonInput.PerPositionVariants[ String( index ) ][ i ];
                let annotation = jsonInput.PerPositionAnnotations[ String( index ) ][ i ];
                samplesData.push( [ j, i, dnaCharToInt( dnaChar ) ] );
                samplesAnnotation.push( [ j, i, annotation ] );
            }
            
            snvCounts[ "Ref." ].push( parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 0 ] ) );
            snvCounts[ "Ref. discarded" ].push( 
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 1 ] ) + 
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 2 ] )
            );
            snvCounts[ "Del." ].push( parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 3 ] ) );
            snvCounts[ "A" ].push( parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 4 ] ) );
            snvCounts[ "A discarded" ].push(
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 5 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 6 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 7 ] )
            );
            snvCounts[ "C" ].push( parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 8 ] ) );
            snvCounts[ "C discarded" ].push(
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 9 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 10 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 11 ] )
            );
            snvCounts[ "G" ].push( parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 12 ] ) );
            snvCounts[ "G discarded" ].push(
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 13 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 14 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 15 ] )
            );
            snvCounts[ "T" ].push( parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 16 ] ) );
            snvCounts[ "T discarded" ].push(
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 17 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 18 ] ) +
                parseInt( jsonInput.PerPositionCounts[ String( index ) ][ 19 ] )
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
                proteinSequenceContent = aaIntToChar( proteinData[ dataIndex ][ 2 ] );
                if ( seriesIndex == 1 ) {
                    // CASE: Item is from samples chart.
                    let index = params.data[ 1 ];
                    sampleName = samplesLabels[ index ];
                    sampleContent = dnaIntToChar( params.data[ 2 ] );
                    let annotationString = samplesAnnotation[ params.dataIndex ][ 2 ];
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
        grid: [
            {   // Grid to display protein sequence information.
                top: '10',
                left: '17%',
                height: '5%',
                width: '77%'
            },
            {   // Grid to display reference sequence information.
                top: '50',
                left: '17%',
                height: '5%',
                width: '77%'
            },
            {   // Grid to display sample content and annotation information.
                top: '90', //'21.5%',
                left: '17%',
                height: '54%',
                width: '77%',
                show: true,
                backgroundColor: '#EEEEEE'
            },
            {   // Grid to display global composition of variants.
                top: 490,
                bottom: 10,
                left: '17.2%',
                width: '76.9%',
                show: true,
                backgroundColor: '#EEEEEE'
            }
        ],
        xAxis: [
            {   // X axis to display reference sequence information.
                type: 'category',
                data: positions,
                splitArea: {
                    show: false
                },
                show: false,
                gridIndex: 1
            },
            {   // X axis to display sample content and annotation information.
                type: 'category',
                data: positions,
                splitArea: {
                    show: false
                },
                show: true,
                gridIndex: 2
            },
            {   // X axis to display global composition of variants.
                type: 'category',
                data: positions,
                splitArea: {
                    show: false
                },
                show: false,
                gridIndex: 3
            },
            {   // X axis to display protein sequence information.
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
            {   // Y axis to display reference sequence information.
                type: 'category',
                data: referenceLabel,
                splitArea: {
                    show: false
                },
                gridIndex: 1
            },
            {   // Y axis to display sample content and annotation information.
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
            {   // Y axis to display global composition of variants.
                type: 'value',
                min: 0,
                max: samplesLabels.length,
                splitArea: {
                    show: false
                },
                show: false,
                gridIndex: 3
            },
            {  // Y axis to display protein sequence information.
                type: 'category',
                data: proteinLabel,
                splitArea: {
                    show: false
                },
                gridIndex: 0
            }
        ],
        visualMap: [
            {   // Visual map to color nucleotide information.
                type: 'piecewise',
                pieces: [
                    // ACGT$.:;-aqygfhtrzcvd
                    { min: 1, max: 1, color: colorSchemes.nucleotides.A }, // A
                    { min: 2, max: 2, color: colorSchemes.nucleotides.C }, // C
                    { min: 3, max: 3, color: colorSchemes.nucleotides.G }, // G
                    { min: 4, max: 4, color: colorSchemes.nucleotides.T }, // T
                    { min: 5, max: 5, color: colorSchemes.nucleotides.nomatch }, // No match with protein sequence.
                    { min: 6, max: 6, color: colorSchemes.nucleotides.ref }, // Reference
                    { min: 7, max: 7, color: colorSchemes.nucleotides.dref }, // Reference, low quality
                    { min: 8, max: 8, color: colorSchemes.nucleotides.dref }, // Reference, low coverage
                    { min: 9, max: 9, color: colorSchemes.nucleotides.del }, // Deletion
                    { min: 10, max: 10, color: colorSchemes.nucleotides.dA }, // A, low quality
                    { min: 11, max: 11, color: colorSchemes.nucleotides.dA }, // A, low coverage
                    { min: 12, max: 12, color: colorSchemes.nucleotides.dA }, // A, low frequency
                    { min: 13, max: 13, color: colorSchemes.nucleotides.dC }, // C, low quality
                    { min: 14, max: 14, color: colorSchemes.nucleotides.dC }, // C, low coverage
                    { min: 15, max: 15, color: colorSchemes.nucleotides.dC }, // C, low frequency
                    { min: 16, max: 16, color: colorSchemes.nucleotides.dG }, // G, low quality
                    { min: 17, max: 17, color: colorSchemes.nucleotides.dG }, // G, low coverage
                    { min: 18, max: 18, color: colorSchemes.nucleotides.dG }, // G, low frequency
                    { min: 19, max: 19, color: colorSchemes.nucleotides.dT }, // T, low quality
                    { min: 20, max: 20, color: colorSchemes.nucleotides.dT }, // T, low coverage
                    { min: 21, max: 21, color: colorSchemes.nucleotides.dT } // T, low frequency
                ],
                seriesIndex: [ 0, 1 ],
                show: false
            },
            {   // Visual map to color amino-acid information.
                type: 'piecewise',
                pieces: [
                    { min: 1, max: 3, color: colorSchemes.aminoacids.H },
                    { min: 4, max: 5, color: colorSchemes.aminoacids.D },
                    { min: 6, max: 10, color: colorSchemes.aminoacids.S },
                    { min: 11, max: 13, color: colorSchemes.aminoacids.F },
                    { min: 14, max: 20, color: colorSchemes.aminoacids.A },
                    { min: 21, max: 21, color: colorSchemes.aminoacids.del },
                    { min: 22, max: 23, color: colorSchemes.aminoacids.nomatch }
                ],
                seriesIndex: [ 13 ],
                show: false
            }
        ],
        legend: { // Disable default legend.
            show: false
        },
        dataZoom: [
            {   // Data zoom to zoom into position intervals.
                type: 'slider',
                xAxisIndex: [ 0, 1, 3 ],
                minValueSpan: 30,
                maxValueSpan: 2900,
                realtime: false,
                throttle: 100,
                top: 482,
                bottom: 14,
                left: '17%',
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
            {   // Data zoom to zoom into single samples.
                type: 'slider',
                yAxisIndex: [ 1 ],
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
                    color: colorSchemes.nucleotides.A
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
                    color: colorSchemes.nucleotides.dA
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
                    color: colorSchemes.nucleotides.C
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
                    color: colorSchemes.nucleotides.dC
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
                    color: colorSchemes.nucleotides.G
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
                    color: colorSchemes.nucleotides.dG
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
                    color: colorSchemes.nucleotides.T
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
                    color: colorSchemes.nucleotides.dT
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
                    color: colorSchemes.nucleotides.del
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
                    color: colorSchemes.nucleotides.dref
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
                    color: colorSchemes.nucleotides.ref
                },
                data: snvCounts[ "Ref." ],
                xAxisIndex: 2,
                yAxisIndex: 2
            },
            {
                name: 'Protein',
                type: 'heatmap',
                data: proteinData,
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
        let proteinSequenceSubstring = proteinSequence.substring( 0, dataIndex + 1 );
        if ( proteinSequenceSubstring.endsWith( "?" ) ) {
            return;
        } else {
            let noMatchCount = ( proteinSequenceSubstring.match( /\?/g ) || [ ] ).length;
            let position = Math.ceil( ( dataIndex + 1 - noMatchCount ) / 3 );
            viewer.removeAllLabels();
            viewer.setStyle({}, {cartoon: {colorfunc: (atom) => { return colorSchemes.mol3d[ atom.resn ] } }});
            viewer.addStyle({resi: [proteinResidueNumberMap[position]]},{stick:{radius:1,color:"#FFBE0B"}});
            viewer.addResLabels({resi: [proteinResidueNumberMap[position]]});
            viewer.render();
        }
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
	let residueIndex = 1;
    let residueNumber = false;
    let oldResidueNumber = false;
    for ( let entry of jsonInput.ProteinStructure.split( "\r\n" ) ) {
        entry = entry.split( " " ).filter( ( str ) => { return /\S/.test( str ) } );
        if ( entry[ 0 ] == "ATOM" ) {
            let resn = parseInt( entry[ 5] );
            if ( ! residueNumber && ! oldResidueNumber ) {
                residueNumber = resn;
                oldResidueNumber = resn;
                proteinResidueNumberMap[ residueIndex ] = residueNumber;
                residueIndex += 1;
            } else {
                residueNumber = resn;
                if ( residueNumber != oldResidueNumber ) {
                    oldResidueNumber = resn;
                    proteinResidueNumberMap[ residueIndex ] = residueNumber;
                    residueIndex += 1;
                }
            }
        }
    }
    let element = $('#proteinView');
    let config = { backgroundColor: '#FFFFFF' };
    viewer = $3Dmol.createViewer( element, config );
    viewer.addModel( jsonInput.ProteinStructure, "pdb" );
    viewer.setStyle({}, {cartoon: {colorfunc: (atom) => { return colorSchemes.mol3d[ atom.resn ] } }});
    viewer.setViewStyle({style:"outline"});
    viewer.zoomTo();
    viewer.render();
    viewer.zoom(1.5, 1000);
}
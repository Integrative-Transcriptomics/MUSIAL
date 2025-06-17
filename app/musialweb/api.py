from musialweb import app
from flask import request, session, send_file
from flask_session import Session
from datetime import timedelta
from dotenv import load_dotenv
from io import StringIO
from cachelib import FileSystemCache
from copy import deepcopy
from sklearn.manifold import trustworthiness, TSNE
from sklearn.cluster import AgglomerativeClustering
import json, zlib, os, subprocess, shutil, random, string, re, base64, copy, traceback, re, time
import igraph as ig
import pandas as pd
import numpy as np
import scipy as sc

""" The working directory inside the Docker container. """
WRK_DIR = "/app/musialweb/"
#WRK_DIR = "./app/musialweb/"
""" Path to the MUSIAL jar executable inside the Docker container. """
MUSIAL = "/app/musialweb/MUSIAL-v2.3.10.jar"
#MUSIAL = "./app/musialweb/MUSIAL-v2.3.10.jar"
load_dotenv() # Load env. variables from local file.
# Set session configuration parameters.
app.config["SECRET_KEY"] = os.getenv("KEY")
app.config["SESSION_COOKIE_NAME"] = "MUSIALweb"
app.config["SESSION_TYPE"] = "cachelib"
app.config["SESSION_PERMANENT"] = True
app.config["PERMANENT_SESSION_LIFETIME"] = timedelta(days=1)
app.config["SESSION_USE_SIGNER"] = True
app.config['SESSION_CACHELIB'] = FileSystemCache(cache_dir = WRK_DIR + 'session/', threshold=0)
app.config["MAX_CONTENT_LENGTH"] = 100 * 1024 * 1024  # Limit content lengths to 100 MB.
"""
Definition of constant session keys:
- DATA: The top-level key for all session data.
- REFERENCE_SEQUENCE: The key for the reference sequence.
- STORAGE: The key for the MUSIAL storage data.
- RECORD_SAMPLES: The key for the sample records; i.e., table columns and records.
- RECORD_FEATURES: The key for the feature records; i.e., table columns and records.
- RECORD_NUCLEOTIDE_VARIANTS: The key for the nucleotide variant records; i.e., table columns and records.
- RECORD_AMINOACID_VARIANTS: The key for the amino acid variant records; i.e., table columns and records.
- META_COMPONENTS: The key for meta-components; i.e., additional data computed from the MUSIAL storage data.
"""
DATA = "data"
REFERENCE_SEQUENCE = "reference_sequence"
STORAGE = "storage"
TABLE_RECORDS = "table_records"
SAMPLES = "samples"
FEATURES = "features"
NUCLEOTIDE_VARIANTS = "nucleotide_variants"
AMINOACID_VARIANTS = "aminoacid_variants"
META_COMPONENTS = "meta_components"
"""
Def. constant misc. variables:
- DEBUG: A boolean indicating whether the server API is in debug mode.
"""
DEBUG = bool(int(os.getenv("DEBUG")))

Session(app) # Start the session.

@app.route("/session/start", methods=["POST"])
def session_start():
    # Clear any previous session data.
    session.clear()
    # Generate unique hex string to use as directory name in the local file system.
    TMP_DIR = _get_tmp_dir_name()
    # Variables to store the output of MUSIAL run.
    stdout = ""
    stderr = ""
    try:
        # Generate directory to store data temporary in the local file system.
        os.makedirs(WRK_DIR + TMP_DIR)
        # Inflate the request data and transform into python dictionary.
        inflated_request_data = zlib.decompress(request.data)
        json_string_request_data = inflated_request_data.decode("utf8")
        json_request_data = json.loads(json_string_request_data)
        # Write reference .fasta to local file and set path in MUSIAL build configuration object.
        session[DATA] = {}
        session[DATA][REFERENCE_SEQUENCE] = json_request_data["referenceSequence"]
        with open(
            WRK_DIR + TMP_DIR + "/reference.fasta", "w+"
        ) as reference_fasta_file:
            reference_fasta_file.write(json_request_data["referenceSequence"])
            del json_request_data["referenceSequence"]
            json_request_data["referenceSequenceFile"] = (
                WRK_DIR + TMP_DIR + "/reference.fasta"
            )
        # Write reference .gff3 to local file and set path in run specification.
        with open(
            WRK_DIR + TMP_DIR + "/reference.gff3", "w+"
        ) as reference_gff3_file:
            reference_gff3_file.write(json_request_data["referenceFeatures"])
            del json_request_data["referenceFeatures"]
            json_request_data["referenceFeaturesFile"] = (
                WRK_DIR + TMP_DIR + "/reference.gff3"
            )
        # For each specified sample, write .vcf to local file and set path in run specification.
        for sample in json_request_data["samples"].keys():
            with open(
                WRK_DIR + TMP_DIR + "/" + sample + ".vcf", "w+"
            ) as sample_vcf_file:
                sample_vcf_file.write(json_request_data["samples"][sample]["vcfFile"])
                json_request_data["samples"][sample]["vcfFile"] = (
                    WRK_DIR + TMP_DIR + "/" + sample + ".vcf"
                )
        # Write the adjusted request (i.e. used as MUSIAL build configuration) to local file.
        with open(
            WRK_DIR + TMP_DIR + "/configuration.json", "w+"
        ) as build_configuration_file:
            json_request_data["output"] = (
                WRK_DIR + TMP_DIR + "/storage.json"
            )
            json.dump(json_request_data, build_configuration_file)
        # Run MUSIAL on the specified data.
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-Xms1G",
                "-Xmx8G",
                "-jar",
                MUSIAL,
                "build",
                "-C",
                WRK_DIR + TMP_DIR + "/configuration.json",
                "-u",
                "-w",
                WRK_DIR + TMP_DIR + "/snpeff/"
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode(encoding="utf-8")
        stderr = stderr.decode(encoding="utf-8")
        # Check, if any error was raised during the MUSIAL run. WARN and INFO entries from BioJava have to be filtered!
        biojava_info_tag = "INFO  org.biojava.nbio"
        biojava_warn_tag = "WARN  org.biojava.nbio"
        if (
            any(
                biojava_info_tag not in e and biojava_warn_tag not in e
                for e in filter(None, stderr.split("\n"))
            )
        ):
            print( stdout )
            print( stderr )
            raise Exception(_remove_ansi(stderr))
        else:
            with open(
                WRK_DIR + TMP_DIR + "/storage.json", "r"
            ) as run_out_file:
                session[DATA][STORAGE] = json.load(run_out_file)
            _process_storage(TMP_DIR)
    # If any error is thrown by the server, set response code to 1 (failed).
    except Exception as e:
        return _format_exception(e), 500
    finally:
        # Remove temporary directory and return response.
        if not DEBUG :
            shutil.rmtree(WRK_DIR + TMP_DIR)
        return "Ok", 200

@app.route("/session/example", methods=["GET"])
def session_example():
    try :
        session.clear()
        with open( WRK_DIR + "/static/resources/exampleSession.zlib", "r" ) as example_session_data :
            json_string_session_data = zlib.decompress( base64.b64decode( example_session_data.read( ) ) ).decode( )
            session[DATA] = json.loads(json_string_session_data)
        with open( WRK_DIR + "/static/resources/H37Rv.fasta", "r" ) as reference_sequence_file :
            session[DATA][REFERENCE_SEQUENCE] = reference_sequence_file.read( )
        
        _process_storage( "" ) # TODO: Remove.
        
        return "Ok", 200
    except Exception as e :
        return _format_exception(e), 500

@app.route("/session/reload", methods=["POST"])
def session_reload():
    try :
        session.clear()
        json_string_session_data = zlib.decompress( base64.b64decode( request.data ) ).decode( )
        session[DATA] = json.loads(json_string_session_data)

        _process_storage( "" ) # TODO: Remove.

        return "Ok", 200
    except Exception as e :
        return _format_exception(e), 500

@app.route("/session/data", methods=["GET"])
def session_data():
    if _validate_session( ):
        _session_data = copy.deepcopy(session[DATA])
        del _session_data[REFERENCE_SEQUENCE]
        return json.dumps(_session_data).replace( "NaN", "null" ), 200
    else :
        return "The requested resource is not available. Please check whether a valid session has been established.", 404

@app.route("/session/save", methods=["GET"])
def session_save():
    try :
        #zlib_compress = zlib.compressobj( 6, zlib.DEFLATED, zlib.MAX_WBITS )
        #compressed_session_bytes = zlib_compress.compress( bytes(json.dumps(session[DATA]), "utf-8") ) + zlib_compress.flush( )
        #encoded_session = base64.b64encode( compressed_session_bytes ).decode("ascii")
        #return encoded_session, 200

        foo = copy.deepcopy(session[DATA])
        del foo["table_records"]
        zlib_compress = zlib.compressobj( 6, zlib.DEFLATED, zlib.MAX_WBITS )
        compressed_session_bytes = zlib_compress.compress( bytes(json.dumps(foo), "utf-8") ) + zlib_compress.flush( )
        encoded_session = base64.b64encode( compressed_session_bytes ).decode("ascii")
        return encoded_session, 200
    except Exception as e :
        return _format_exception(e), 500

@app.route("/compute/correlation", methods=["POST"])
def compute_correlation():
    try :
        inflated_request_data = zlib.decompress(request.data)
        json_string_request_data = inflated_request_data.decode("utf8")
        json_request_data = json.loads(json_string_request_data)
        col1 = json_request_data["col1"]
        col2 = json_request_data["col2"]
        test_type = json_request_data["test_type"]
        entry_type = json_request_data["entry_type"]
        df = pd.DataFrame.from_records( session[DATA]["record_" + entry_type]["records"] )
        if test_type == "pearsonr":
            t, p = sc.stats.pearsonr(df[col1].to_numpy(), df[col2].to_numpy())
        elif test_type == "spearmanr":
            t, p = sc.stats.spearmanr(df[col1].to_numpy(), df[col2].to_numpy())
        elif test_type == "kendalltau":
            t, p = sc.stats.kendalltau(df[col1].to_numpy(), df[col2].to_numpy())
        elif test_type == "cramer":
            t = sc.stats.contingency.association(
                df.groupby([col1, col2])
                .size()
                .unstack()
                .replace(np.nan, 0)
                .astype(int)
                .to_numpy(),
                method="cramer",
            )
            p = None
        else:
            raise Exception("Test type not implemented.")
        return json.dumps({"t": round( t, 4 ) if t != None else t, "p": round( p, 4 ) if p != None else p }).replace( "NaN", "null" ), 200
    except Exception as e:
        return _format_exception(e), 500

@app.route("/download/sequences", methods=["POST"])
def download_sequences():
    # Generate unique hex string to use as directory name in the local file system.
    TMP_DIR = _get_tmp_dir_name()
    try:
        # Inflate the request data and transform into python dictionary.
        inflated_request_data = zlib.decompress(request.data)
        json_string_request_data = inflated_request_data.decode("utf8")
        json_request_data = json.loads(json_string_request_data)
        # Generate directory to store data temporary in the local file system.
        os.makedirs( WRK_DIR + TMP_DIR )
        with open( WRK_DIR + TMP_DIR + "/storage.json", "w+" ) as session_storage:
            session_storage.write(json.dumps(session[DATA][STORAGE]))
        if json_request_data[ "content" ] == "nucleotide" and json_request_data[ "conserved" ] :
            with open( WRK_DIR + TMP_DIR + "/reference.fasta", "w+" ) as nucleotide_reference:
                nucleotide_reference.write(session[DATA][REFERENCE_SEQUENCE])
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-Xms1G",
                "-Xmx8G",
                "-jar",
                MUSIAL,
                "export_sequence",
                "-a" if json_request_data[ "align" ] else "",
                "-c",
                json_request_data[ "content" ],
                "-F",
                json_request_data[ "feature" ],
                "-m" if json_request_data[ "merge" ] else "",
                "-I",
                WRK_DIR + TMP_DIR + "/storage.json",
                "-k" if json_request_data[ "conserved" ] else "",
                "-O",
                WRK_DIR + TMP_DIR + "/out.fasta",
                "-r" if ( json_request_data[ "conserved" ] and json_request_data[ "content" ] == "nucleotide" ) else "",
                WRK_DIR + TMP_DIR + "/reference.fasta" if ( json_request_data[ "conserved" ] and json_request_data[ "content" ] == "nucleotide" ) else "",
                "-s",
                *json_request_data[ "samples" ]
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode(encoding="utf-8")
        stderr = stderr.decode(encoding="utf-8")
        if stderr != "":
            raise Exception(_remove_ansi(stdout) + "\n" + _remove_ansi(stderr))
        else :
            return send_file(
                TMP_DIR + "/out.fasta",
                as_attachment=True,
            ), 200
    except Exception as e:
        return _format_exception(e), 500
    finally:
        if not DEBUG :
            shutil.rmtree(TMP_DIR)

@app.route("/download/variants", methods=["POST"])
def download_variants():
    # Generate unique hex string to use as directory name in the local file system.
    TMP_DIR = _get_tmp_dir_name()
    try:
        # Inflate the request data and transform into python dictionary.
        inflated_request_data = zlib.decompress(request.data)
        json_string_request_data = inflated_request_data.decode("utf8")
        json_request_data = json.loads(json_string_request_data)
        # Generate directory to store data temporary in the local file system.
        os.makedirs( WRK_DIR + TMP_DIR )
        with open( WRK_DIR + TMP_DIR + "/storage.json", "w+" ) as session_storage:
            session_storage.write(json.dumps(session[DATA][STORAGE]))
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-Xms1G",
                "-Xmx8G",
                "-jar",
                MUSIAL,
                "export_table",
                "-c",
                json_request_data[ "content" ],
                "-F",
                json_request_data[ "feature" ],
                "-I",
                WRK_DIR + TMP_DIR + "/storage.json",
                "-O",
                WRK_DIR + TMP_DIR + "/out.tsv",
                "-s",
                *json_request_data[ "samples" ]
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode(encoding="utf-8")
        stderr = stderr.decode(encoding="utf-8")
        if stderr != "":
            raise Exception(_remove_ansi(stdout) + "\n" + _remove_ansi(stderr))
        else :
            return send_file(
                TMP_DIR + "/out.tsv",
                as_attachment=True,
            ), 200
    except Exception as e:
        return _format_exception(e), 500
    finally:
        if not DEBUG :
            shutil.rmtree(TMP_DIR)

@app.route("/download/sequence_typing", methods=["GET"])
def download_sequence_typing():
    # TODO: Implement this in MUSIAL.
    TMP_DIR = _get_tmp_dir_name()
    try :
        sequence_type_idx = 1
        sequence_types = { }
        content = "sample\t" + \
            "\t".join( [ "allele_" + str( _ ) for _ in session[DATA][STORAGE]["features"].keys( ) ] ) + \
            "\tsequence_type\n"
        for sample in session[DATA][STORAGE]["samples"].keys( ) :
            content += sample + "\t"
            sample_allele_profile = tuple(
                [ "0" if _ == "reference" else _ for _ in 
                    [ _.strip( "A" ) for _ in session[DATA][STORAGE]["samples"][sample]["allele"].values( ) ]
                ]
            )
            if not sample_allele_profile in sequence_types :
                sequence_types[ sample_allele_profile ] = sequence_type_idx
                sequence_type_idx += 1
            content += "\t".join( sample_allele_profile ) + "\t" + str( sequence_types[ sample_allele_profile ] ) + "\n"
        os.makedirs( WRK_DIR + TMP_DIR )
        with open( WRK_DIR + TMP_DIR + "/sequence_typing.tsv", "w+" ) as sequence_typing_file :
            sequence_typing_file.write( content )
        return send_file(
                TMP_DIR + "/sequence_typing.tsv",
                as_attachment=True,
            ), 200
    except Exception as e:
        return _format_exception(e), 500
    finally:
        if not DEBUG :
            shutil.rmtree(TMP_DIR)

def _validate_session():
    """
    Validates the current session by checking whether all required data is available.

    Returns:
        bool: True if the session is valid; otherwise, False.
    """
    return DATA in session and REFERENCE_SEQUENCE in session[DATA] and STORAGE in session[DATA] and TABLE_RECORDS in session[DATA] and all( _ in session[DATA][TABLE_RECORDS] for _ in (SAMPLES, FEATURES, NUCLEOTIDE_VARIANTS, AMINOACID_VARIANTS) )
    
def _process_storage(storage_dir: str):
    """
    Re-process the storage data by ... TODO

    Args:
        storage_dir (str): The directory path where the storage data is located (without the PATH_PREFIX).

    Raises:
        Exception: If no storage data is available for the session.

    Returns:
        None
    """
    if not STORAGE in session[DATA]:
        raise Exception("No storage data available for session. Has a session been started?")
    
    if not TABLE_RECORDS in session[DATA]:
        session[DATA][TABLE_RECORDS] = { }

    # (i) Run MUSIAL on the specified data to view samples.
    if not SAMPLES in session[DATA][TABLE_RECORDS]:
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-Xms1G",
                "-Xmx8G",
                "-jar",
                MUSIAL,
                "view_samples",
                "-I",
                WRK_DIR + storage_dir + "/storage.json",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode()
        stderr = stderr.decode()
        session[DATA][TABLE_RECORDS][SAMPLES] = _get_records_from_view_output(stdout)

    # (ii) Run MUSIAL on the specified data to view features.
    if not FEATURES in session[DATA][TABLE_RECORDS]:
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-Xms1G",
                "-Xmx8G",
                "-jar",
                MUSIAL,
                "view_features",
                "-I",
                WRK_DIR + storage_dir + "/storage.json",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode(encoding="utf-8")
        stderr = stderr.decode(encoding="utf-8")
        session[DATA][TABLE_RECORDS][FEATURES] = _get_records_from_view_output(stdout)

    # (iii) Run MUSIAL on the specified data to view variants.
    if not NUCLEOTIDE_VARIANTS in session[DATA][TABLE_RECORDS]:
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-Xms1G",
                "-Xmx8G",
                "-jar",
                MUSIAL,
                "view_variants",
                "-I",
                WRK_DIR + storage_dir + "/storage.json",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode(encoding="utf-8")
        stderr = stderr.decode(encoding="utf-8")
        session[DATA][TABLE_RECORDS][NUCLEOTIDE_VARIANTS] = _get_records_from_view_output(stdout)
    if not AMINOACID_VARIANTS in session[DATA][TABLE_RECORDS]:
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-Xms1G",
                "-Xmx8G",
                "-jar",
                MUSIAL,
                "view_variants",
                "-c",
                "aminoacid",
                "-I",
                WRK_DIR + storage_dir + "/storage.json",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode(encoding="utf-8")
        stderr = stderr.decode(encoding="utf-8")
        session[DATA][TABLE_RECORDS][AMINOACID_VARIANTS] = _get_records_from_view_output(stdout)

    # (iv) Add meta-components to session data.
    if not META_COMPONENTS in session[DATA]:
        session[DATA][META_COMPONENTS] = {}

    # (iv.1) Run variant counts clustering:
    #if not "variant_counts_clustering" in session[DATA][META_COMPONENTS]:
    session[DATA][META_COMPONENTS]["variant_counts_clustering"] = _compute_variant_counts_clustering()

    # (iv.2) Compute sample distance clustering:
    #if not "sample_distance_manifold" in session[DATA][META_COMPONENTS]:
    session[DATA][META_COMPONENTS]["sample_distance_manifold"] = _compute_sample_distance_manifold()

    # (iv.3) Prepare empty meta-component for gene networks:
    #if not "allele_networks" in session[DATA][META_COMPONENTS]:
    session[DATA][META_COMPONENTS]["allele_networks"] = _compute_allele_networks()

def _compute_variant_counts_clustering() :
    _start = time.time( )
    S = session[DATA][STORAGE]
    variant_counts = { }
    for entry in S["samples"].values() :
        sample_name = entry[ "name" ]
        variant_counts.setdefault( sample_name, { } )
        for feature_name in S["features"].keys() :
            variant_counts[ sample_name ].setdefault( feature_name, 0 )
            if feature_name in entry[ "calls" ] :
                for call in entry[ "calls" ][ feature_name ] :
                    callinfo = call.split( "," )[0]
                    if callinfo != "0" and callinfo != "?" :
                        variant_counts[sample_name][feature_name] += 1
    variant_counts_matrix = [ [ sn, fn, c ] for sn in variant_counts.keys( ) for fn, c in variant_counts[ sn ].items( ) ]
    variant_counts_df = pd.DataFrame( variant_counts_matrix ).pivot( index = 0, columns = 1, values = 2 )
    variant_counts_df.index.name = None
    variant_counts_df.columns.name = None
    variant_counts_df.fillna( 0, inplace = True )

    _sample_names = list( variant_counts_df.index )
    _sample_hierachical_clustering = AgglomerativeClustering( n_clusters = int( np.ceil( np.log2( len( _sample_names ) ) ) ) )
    if len( _sample_names ) > 1 :
        _sample_hierachical_clustering.fit( variant_counts_df.to_numpy( ) )
    else :
        _sample_hierachical_clustering.labels_ = [ 0 ]

    _feature_names = list( variant_counts_df.columns )
    _feature_hierachical_clustering = AgglomerativeClustering( n_clusters = int( np.ceil( len( _feature_names ) / 2 ) ) )
    if len( _feature_names) > 1 :
        _feature_hierachical_clustering.fit( np.transpose( variant_counts_df.to_numpy( ) ) )
    else :
        _feature_hierachical_clustering.labels_ = [ 0 ]

    _end = time.time( )
    _debug( "Runtime _compute_variant_counts_clustering: " + str( _end - _start ) )
    return {
        "description": "(Independent) clustering of samples (across features) and features (across samples) by variant counts.",
        "counts": variant_counts_matrix,
        "sample_labels": { k: int( v ) for k, v in zip( _sample_names, _sample_hierachical_clustering.labels_ ) },
        "feature_labels": { k: int( v ) for k, v in zip( _feature_names, _feature_hierachical_clustering.labels_ ) }
    }

def _compute_sample_distance_manifold() :
    _start = time.time( )
    S = session[DATA][STORAGE]
    # 1. Infer allele profile for samples:
    allele_profiles = { }
    # 1.1 Add reference profile first by default.
    reference_profile = [ ]
    for feature in S["features"].keys() :
        reference_profile.append( feature + ".reference" )
    allele_profiles.setdefault( "/".join( reference_profile ), {
        "samples": [ ],
    } )
    # 1.2 Infer allele profile for each sample.
    for name in S["samples"].keys( ) :
        alleles = [ ]
        for feature, allele in S["samples"][name]["allele"].items( ) :
            alleles.append( feature + "." + allele )
        allele_profile = "/".join( alleles )
        if not allele_profile in allele_profiles :
            allele_profiles.setdefault( allele_profile, {
                "samples": [ ],
            } )
        allele_profiles[ allele_profile ][ "samples" ].append( name )

    # 2. Compute allele profile distances:
    distance_matrix = [ ]
    allele_profile_distances = { }
    def allele_profile_distance( ap1, ap2 ) :
        score = 0
        for fa1, fa2 in zip( ap1.split( "/" ), ap2.split( "/" ) ) :
            if fa1 != fa2 :
                feature = fa1.split( "." )[0]
                a1 = fa1.split( "." )[1]
                a2 = fa2.split( "." )[1]
                if a1 == "reference" and a2 == "reference" :
                    continue
                else :
                    a1Variants = set( [ ] ) if a1 == "reference" else set( S[ "features" ][ feature ][ "alleles" ][ a1 ][ "variants" ].split( "," ) )
                    a2Variants = set( [ ] ) if a2 == "reference" else set( S[ "features" ][ feature ][ "alleles" ][ a2 ][ "variants" ].split( "," ) )
                    symmetric_difference = a1Variants.symmetric_difference( a2Variants )
                    for variant in symmetric_difference :
                        pos, alt = variant.split( ":" )
                        impact = S[ "features" ][ feature ][ "nucleotideVariants" ][ pos ][ alt ][ "info" ].get( "snpeff_Impact", False )
                        if impact :
                            if impact == "HIGH" :
                                score += 1
                            elif impact == "MODERATE" :
                                score += 1
                            elif impact == "LOW" :
                                score += 1
                            elif impact == "MODIFIER" :
                                score += 1
        return score
    for ap1 in allele_profiles.keys( ) :
        _ = [ ]
        for ap2 in allele_profiles.keys( ) :
            if not ( ap1, ap2 ) in allele_profile_distances and not ( ap2, ap1 ) in allele_profile_distances :
                if ap1 == ap2 :
                    allele_profile_distances[ ( ap1, ap2 ) ] = 0
                else :
                    allele_profile_distances[ ( ap1, ap2 ) ] = allele_profile_distances[ ( ap2, ap1 ) ] = allele_profile_distance( ap1, ap2 )
            _.append( allele_profile_distances[ ( ap1, ap2 ) ] )
        distance_matrix.append( _ )
    distance_matrix = np.array( distance_matrix )

    # 3. Run dimensionality reduction on distance matrix:
    allele_profile_labels = list( allele_profiles.keys( ) )
    N = len( allele_profile_labels )
    allele_profiles_manifold = TSNE(
        n_components = 2,
        perplexity = np.ceil( N / 8 ),
        early_exaggeration = 12,
        metric = "precomputed",
        init = "random",
        random_state = 0,
        angle = 0.1
    )
    distance_matrix_projection = allele_profiles_manifold.fit_transform( distance_matrix )

    # 4. Re-expand data to sample level:
    samples = { }
    for apIndex in range( len( allele_profile_labels ) ) :
        for sample in allele_profiles[ allele_profile_labels[ apIndex ] ][ "samples" ] :
            samples.setdefault(
                sample,
                {
                    "projection": distance_matrix_projection[ apIndex ].tolist( ),
                    "apIndex": apIndex
                }
            )
    _end = time.time( )
    _debug( "Runtime _compute_sample_distance_manifold: " + str( _end - _start ) )
    return {
        "description": "MDS (per sample) of pairwise distances between allele profiles of samples.",
        "samples": samples,
        "mean_nnd": sc.spatial.cKDTree( distance_matrix_projection ).query( distance_matrix_projection, k = 2 )[0][ :, 1 ].mean( ),
        "trustworthiness": trustworthiness( distance_matrix, distance_matrix_projection, metric = "precomputed", n_neighbors = int( np.ceil( N / 3 ) ) )
        #"allele_profile_labels": allele_profile_labels,
        #"allele_profiles": allele_profiles,
    }

def _compute_allele_networks() :
    _start = time.time( )
    allele_networks = { }
    _session = session[DATA][STORAGE]
    for feature_name in _session[ "features" ].keys() :
        _feature = _session[ "features" ][feature_name]
        alleles_variants_dict = dict(
            sorted(
                {
                    allele_name : set( [ ] ) if allele_name == "reference" else set( _session[ "features" ][ feature_name ][ "alleles" ][ allele_name ][ "variants" ].split( "," ) )
                    for allele_name in list( _session[ "features" ][ feature_name ][ "alleles" ].keys( ) )
                }.items( ),
                key = lambda _ : len( _[1] ),
                reverse = True
            )
        )
        alleles_linkage = {
            "reference": [ ]
        }
        #alleles_is_sink = dict( zip( alleles_variants_dict.keys( ), [ True ] * len( alleles_variants_dict ) ) )
        #alleles_is_sink[ "reference" ] = False
        for allele_name in alleles_variants_dict.keys( ) :
            if allele_name == "reference" :
                continue
            allele_variants = copy.deepcopy( alleles_variants_dict[ allele_name ] )
            other_allele_names = set( alleles_variants_dict.keys( ) )
            other_allele_names.remove( "reference" )
            other_allele_names.remove( allele_name )
            other_allele_names = list( filter( lambda _ : len( alleles_variants_dict[ _ ] ) < len( allele_variants ), other_allele_names ) )
            iterate = True
            is_linked = False
            while iterate :
                shared_variants = { }
                for other_allele_name in other_allele_names :
                    intersection = allele_variants.intersection( alleles_variants_dict[ other_allele_name ] )
                    if len( intersection ) > 0 :
                        symmetric_difference = alleles_variants_dict[ allele_name ].symmetric_difference( alleles_variants_dict[ other_allele_name ] )
                        shared_variants[ other_allele_name ] = ( symmetric_difference, intersection )
                shared_variants = sorted( shared_variants.items( ), key = lambda _ : len( _[ 1 ][ 0 ] ) - len( _[ 1 ][ 1 ] ) )
                if len( shared_variants ) > 0 :
                    _ = shared_variants[ 0 ]
                    alleles_linkage.setdefault( allele_name, [ ] )
                    alleles_linkage[ allele_name ].append( ( _[ 0 ], len( _[ 1 ][ 0 ] ) ) )
                    #alleles_is_sink[ _[ 0 ] ] = False
                    allele_variants = allele_variants.difference( _[ 1 ][ 1 ] )
                    other_allele_names.remove( _[ 0 ] )
                    is_linked = True
                else :
                    iterate = False
            if not is_linked :
                alleles_linkage.setdefault( allele_name, [ ] )
                alleles_linkage[ allele_name ].append( ( "reference", len( allele_variants ) ) )

        network = ig.Graph( directed = True )
        for allele_name in alleles_linkage.keys( ) :
            network.add_vertex(
                name = allele_name,
                number_of_samples = len( _feature["alleles"][allele_name]["occurrence"] ) if allele_name in _feature["alleles"] else 0,
                proteoform = _feature["alleles"][allele_name]["info"]["proteoform"] if allele_name in _feature["alleles"] and "proteoform" in _feature["alleles"][allele_name]["info"] else "N/A",
                #is_sink = alleles_is_sink[ allele_name ]
            )
        for allele_name, _ in alleles_linkage.items( ) :
            for link in _ :
                network.add_edge(
                    source = link[ 0 ],
                    target = allele_name,
                    source_name = link[ 0 ],
                    target_name = allele_name,
                    weight = link[ 1 ]
                )

        scaledWeights = [ _[ "weight" ] * network.vs[ _.target ][ "number_of_samples" ] for _ in network.es ]
        layout = network.layout_fruchterman_reingold(
            weights = [ 1 / _ for _ in scaledWeights ],
            seed = network.layout_circle( ).coords
        )
        network_dict = { "nodes": [ ], "edges": [ ] }
        for node in network.vs :
            network_dict["nodes"].append( {
                "name": node["name"],
                "value": layout[node.index],
                "_info": [ node["number_of_samples"], node["proteoform"] ]
            } )
        for edge in network.es :
            network_dict["edges"].append( {
                "source": edge.source,
                "target": edge.target,
                "value": edge["weight"],
                "_info": [ edge["source_name"], edge["target_name"] ]
            } )
        allele_networks[ feature_name ] = network_dict

    _end = time.time( )
    _debug( "Runtime _compute_allele_networks: " + str( _end - _start ) )
    return allele_networks

def _get_records_from_view_output(stdout: str) -> dict:
    """
    Extracts records from the stdout of a MUSIAL view command.

    Args:
        stdout (str): The stdout output of the MUSIAL view command.
        count_categories (bool, optional): Whether to count categories. Defaults to False.

    Returns:
        dict: A dictionary containing the columns, records, and counts (if count_categories is True).

    """
    # Remove MUSIAL view specific comment lines from output.
    records = _remove_ansi(stdout).split("\n")[4:-3]
    columns = records[0].split("\t")
    df = pd.read_csv(StringIO("\n".join(records)), sep="\t")
    records = df.to_dict(orient="records")
    return {"columns": columns, "rows": records, "filters": [] }

def _get_tmp_dir_name():
    """
    Generates the name for a temporary working directory.

    Returns:
        str: tmp_X with X being a random string consisting of letters and digits of length 10.
    """
    return 'tmp_' + ''.join(
        random.SystemRandom().choice(string.ascii_letters + string.digits)
        for _ in range(10)
    )

def _remove_ansi(text: str) -> str:
    """
    Removes ANSI escape sequences from the given text.

    Args:
        text (str): The text containing ANSI escape sequences.

    Returns:
        str: The text with ANSI escape sequences removed.
    """
    ansi_remove_expression = re.compile(r"(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]")
    return ansi_remove_expression.sub('', text)
 
def _format_exception(exc: Exception) -> str:
    """
    Formats the given exception into a string.

    If `DEBUG` in the .env file is True/1, the traceback will be printed with color highlighting; otherwise, only the exception message will be formatted.

    Args:
        e (str): The exception to be formatted.

    Returns:
        str: The formatted exception as a string.
    """
    if DEBUG: # Print the exception traceback with color highlighting.
        print( '\u001b[31m[Exception]\u001b[0m\n' )
        print( ''.join( traceback.format_exception(RuntimeError, exc, exc.__traceback__) ) )
    # Format only the exception message.
    return ''.join(traceback.format_exception_only(RuntimeError, exc)).strip()

def _debug(text: str):
    if DEBUG:
        print( '\u001b[44m[DEBUG] ' + text + '\u001b[0m' )
#
# Compute genome MLST using https://github.com/tseemann/mlst
#

use Bio::KBase::AppService::AppConfig qw(application_backend_dir);
use P3DataAPI;
use strict;
use Getopt::Long::Descriptive;  
use File::Temp;
use File::Copy;
use File::Slurp;
use GenomeTypeObject;
use Cwd;
use Data::Dumper;
use Time::HiRes qw(gettimeofday);
use JSON::XS;
use Text::CSV_XS qw(csv);
use IPC::Run qw(run);

my ($opt, $usage) = describe_options("%c %o [< in] [> out]",
                     ["in|i=s", "Input GTO"],
                     ["out|o=s", "Output GTO"],
                     ["parallel=i", "Number of threads to use", { default => 1 }],
                     ["dry_run|n", "Dry run - print commands without executing"],
                     ["help|h", "Print this help message"]);
                     
print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 0);

chomp(my $hostname = `hostname -f`);

# helper sub to run or dry-run a command
sub run_cmd {
    my ($cmd_ref, $label) = @_;
    print STDERR "Invoke $label command: @$cmd_ref\n";
    if ($opt->dry_run) {
        print STDERR "DRY RUN - would execute: @$cmd_ref\n";
        return 1;
    }
    my($stdout, $stderr);
    my $ok = run($cmd_ref,
                 ">", \$stdout,
                 "2>", \$stderr);
    $ok or die "Error $? running $label: @$cmd_ref\n";
    return $ok;
}

my $gto_in;

if ($opt->in)
{
    $gto_in = GenomeTypeObject->new({file => $opt->in});
    $gto_in or die "Error reading GTO from " . $opt->in . "\n";
}
else
{
    $gto_in = GenomeTypeObject->new({file => \*STDIN});
    $gto_in or die "Error reading GTO from standard input\n";
}

my $in_file = File::Temp->new();
$gto_in->write_contigs_to_file($in_file);

my %taxon_to_schema = (
    470     => { schema => "Acinetobacter_baumannii",                              version => "1.0" },
    1392    => { schema => "Bacillus_anthracis",                                   version => "1.0" },
    520     => { schema => "Bordetella_pertussis",                                 version => "1.0" },
    235     => { schema => "Brucella_spp",                                         version => "1.0" },
    36855   => { schema => "Brucella_spp",                                         version => "1.0" },
    120577  => { schema => "Brucella_spp",                                         version => "1.0" },
    1218315 => { schema => "Brucella_spp",                                         version => "1.0" },
    29459   => { schema => "Brucella_melitensis",                                  version => "1.0" },
    444163  => { schema => "Brucella_spp",                                         version => "1.0" },
    29460   => { schema => "Brucella_spp",                                         version => "1.0" },
    236     => { schema => "Brucella_spp",                                         version => "1.0" },
    120576  => { schema => "Brucella_spp",                                         version => "1.0" },
    29461   => { schema => "Brucella_spp",                                         version => "1.0" },
    13373   => { schema => "Burkholderia_mallei_fli",                              version => "1.0" },
    111527  => { schema => "Burkholderia_pseudomallei",                            version => "1.1" },
    197     => { schema => "Campylobacter_jejuni_coli",                            version => "1.3" },
    195     => { schema => "Campylobacter_jejuni_coli",                            version => "1.3" },
    1496    => { schema => "Clostridioides_difficile",                             version => "2.0" },
    1502    => { schema => "Clostridium_perfringens",                              version => "1.0" },
    1717    => { schema => "Corynebacterium_diphtheriae",                          version => "1.0" },
    1719    => { schema => "Corynebacterium_pseudotuberculosis",                   version => "1.0" },
    413503  => { schema => "Cronobacter_sakazakii_malonaticus",                   version => "1.0" },
    28141   => { schema => "Cronobacter_sakazakii_malonaticus",                   version => "1.0" },
    1351    => { schema => "Enterococcus_faecalis",                                version => "1.0" },
    1352    => { schema => "Enterococcus_faecium",                                 version => "1.0" },
    562     => { schema => "Escherichia_coli",                                     version => "1.0" },
    263     => { schema => "Francisella_tularensis",                               version => "1.0" },
    2058152 => { schema => "Klebsiella_oxytoca_grimontii_michiganensis_pasteurii", version => "1.0" },
    1134687 => { schema => "Klebsiella_oxytoca_grimontii_michiganensis_pasteurii", version => "1.0" },
    571     => { schema => "Klebsiella_oxytoca_grimontii_michiganensis_pasteurii", version => "1.0" },
    2587529 => { schema => "Klebsiella_oxytoca_grimontii_michiganensis_pasteurii", version => "1.0" },
    573     => { schema => "Klebsiella_pneumoniae_variicola_quasipneumoniae",      version => "1.0" },
    1463165 => { schema => "Klebsiella_pneumoniae_variicola_quasipneumoniae",      version => "1.0" },
    244366  => { schema => "Klebsiella_pneumoniae_variicola_quasipneumoniae",      version => "1.0" },
    446     => { schema => "Legionella_pneumophila",                               version => "1.0" },
    1639    => { schema => "Listeria_monocytogenes",                               version => "1.0" },
    33894   => { schema => "Mycobacterium_tuberculosis_bovis_africanum_canettii",  version => "2.1" },
    1765    => { schema => "Mycobacterium_tuberculosis_bovis_africanum_canettii",  version => "2.1" },
    78331   => { schema => "Mycobacterium_tuberculosis_bovis_africanum_canettii",  version => "2.1" },
    77643   => { schema => "Mycobacterium_tuberculosis_bovis_africanum_canettii",  version => "2.1" },
    36809   => { schema => "Mycobacteroides_abscessus",                            version => "1.0" },
    2096    => { schema => "Mycoplasma_gallisepticum",                             version => "1.0" },
    1464    => { schema => "Paenibacillus_larvae",                                 version => "1.0" },
    287     => { schema => "Pseudomonas_aeruginosa",                               version => "1.0" },
    28901   => { schema => "Salmonella_enterica",                                  version => "2.0" },
    615     => { schema => "Serratia_marcescens",                                  version => "1.0" },
    1280    => { schema => "Staphylococcus_aureus",                                version => "1.3" },
    29388   => { schema => "Staphylococcus_capitis",                               version => "1.0" },
    1314    => { schema => "Spyogenes",                                             version => "1.0" },
    630     => { schema => "Yersinia_enterocolitica",                              version => "1.0" },
);


# grab all lineage 
my $api = P3DataAPI->new();
my $taxon_id = $gto_in->{ncbi_taxonomy_id};
my @res = $api->query("taxonomy", ['eq', 'taxon_id', $taxon_id], ['select', 'taxon_name', 'lineage_ids', 'lineage_names']);

# Query schema_map hash (taxon_id => schema dir name)
my ($dir_name, $schema_version);
foreach my $lineage_id (reverse @{$res[0]->{lineage_ids} // [] }) {
    $lineage_id = int($lineage_id);
    if (exists $taxon_to_schema{$lineage_id}) {
        $dir_name       = $taxon_to_schema{$lineage_id}{schema};
        $schema_version = $taxon_to_schema{$lineage_id}{version};
        last;
    }
}

print STDERR "dir_name: $dir_name\n";

# based on the schema name, determine there is a cgmlst scheme to use
if ($dir_name) {
    my $tmp_dir = File::Temp->newdir(CLEANUP => 0);
    my $clean_fasta_dir = "$tmp_dir/clean_fastas";
    my $allele_call_out = "$tmp_dir/new_genomes_allele_call/";
    my $schema_path = application_backend_dir . "/CoreGenomeMLST/chewbbaca_schemas/" . $dir_name;
    
    mkdir $clean_fasta_dir or die "Cannot create clean fastas directory: $!\n";
    
    # rename to be certain it ends with 'fasta' for chewBBACA
    run_cmd(["cp", "$in_file", "$clean_fasta_dir/input.fasta"], "Add extension");

    # new genome allele call
    my @allele_call_cmd = (
        "chewBBACA.py", "AlleleCall",
        "--input-files", $clean_fasta_dir,
        "--schema-directory", $schema_path,
        "--output-directory", $allele_call_out,
        "--cpu", "1",
        "--output-unclassified",
        "--output-missing",
        "--output-novel",
        "--no-inferred"
    );
    run_cmd(\@allele_call_cmd, "Allele Call");
    
    # Check percent of loci have an exact allele match
    my $cluster_threshold  = 85;  # run clustering
    my $qc_threshold     = 70;  # good vs poor
    print STDERR "Thresholds: qc_threshold=${qc_threshold}% (good/poor), cluster_threshold=${cluster_threshold}% (run clustering)\n";
    my $new_allele_call = "$allele_call_out/results_alleles.tsv";
    my ($total, $exact, $pct) = (0, 0, 0);
    my $qc = "poor";              # good|poor
    my $do_cluster = 0;           # 1 => run clustering

    if ($opt->dry_run) {
        print STDERR "DRY RUN - would check allele call quality from $new_allele_call\n";
        $total = 100;
        $exact = 100;
        $pct = 100; 
    } else {
        open(my $check_fh, '<', $new_allele_call)
            or die "Cannot open $new_allele_call: $!";
        my $hdr = <$check_fh>;   # skip header
        my $data = <$check_fh>;  # data line
        close($check_fh);

        chomp $data;
        my @vals = split(/\t/, $data);
        shift @vals;  # remove first column (filename)

        $total = scalar @vals;
        print STDERR "Total loci in allele call: $total\n";
        #$exact = grep { /^\d+$/ } @vals;
        $exact = grep { my $v = $_; $v =~ s/^\s+|\s+$//g; $v =~ /^\d+$/ || $v =~ /^INF-\d+$/ } @vals;
        print STDERR "Exact allele matches (integers or INF-###): $exact  |  Non-exact (LNF/missing/etc): " . ($total - $exact) . "\n";
        $pct   = $total > 0 ? ($exact / $total) * 100 : 0;
        }

        $qc = ($pct >= $qc_threshold) ? "good" : "poor";
        $do_cluster = ($pct >= $cluster_threshold) ? 1 : 0;

        print STDERR sprintf(
        "Allele call quality: %d / %d exact matches (%.1f%%) => qc=%s, cluster=%s\n",
        $exact, $total, $pct, $qc, ($do_cluster ? "yes" : "no")
        );
    
        # Parse allele calls into a comma string ONCE (always saved)
        my $allele_call_string;
        if ($opt->dry_run) {
            print STDERR "DRY RUN - would parse allele calls from $new_allele_call\n";
            $allele_call_string = "DRY_RUN_ALLELE_CALLS";
        } else {
            open(my $allele_fh, '<', $new_allele_call) or die "Cannot open $new_allele_call: $!";
            my $header_line = <$allele_fh>;
            my $data_line   = <$allele_fh>;
            close($allele_fh);

            die "No allele-call data line found when parsing string from $new_allele_call\n"
                unless defined($data_line);

            chomp $data_line;
            $data_line =~ s/\t/,/g;        # tabs -> commas
            $data_line =~ s/^[^,]*,//;     # remove first column (filename)
            $allele_call_string = $data_line;
            print STDERR "Allele call string: $allele_call_string\n";
        }

         # Build analysis event once
        my $event = {
            tool_name      => "p3x-compute-cgmlst",
            parameters     => [map { "$_" } @allele_call_cmd],
            execution_time => scalar gettimeofday,
            hostname       => $hostname,
        };
        my $event_id = $gto_in->add_analysis_event($event);

        my $sequence_typing = {
            allele_profile  => $allele_call_string,
            scheme_name   => $dir_name,
            scheme_version => $schema_version,
            event_id      => $event_id,
            loci_total    => $total,
            loci_called   => $exact,
            loci_missing  => $total - $exact,
            pct_called    => $pct,
            qc            => $qc, 
        };

            if (!$do_cluster) {
        print STDERR "Skipping clustering — allele call quality below ${cluster_threshold}%.\n";
        push(@{$gto_in->{sequence_types}}, $sequence_typing);
        } else {

            # Join new allele call with master allele call using chewBBACA JoinProfiles
            my $master_table = application_backend_dir . "/CoreGenomeMLST/precomputed_clusters/refs/" . lc($dir_name) . "_11_25_2025_joined.tsv";
            my $master_joined = "$tmp_dir/master_joined.tsv";
            my $copy_of_master = "$tmp_dir/master_copy.tsv";

            run_cmd(["cp",  $master_table, $copy_of_master], "Copy Master Table");

            run_cmd([
                "chewBBACA.py", "JoinProfiles",
                "--profiles", $copy_of_master, $new_allele_call,
                "--output-file", $master_joined
            ], "Join Profiles");

            run_cmd([
                "core-genome-mlst-utils",
                "clean-allelic-profile", $master_joined
            ], "Clean Allele Call");

            my $precomputed_clusters_dir  = application_backend_dir . "/CoreGenomeMLST/precomputed_clusters/";
            my $precomputed_clusters_path_npz = $precomputed_clusters_dir . lc($dir_name) . ".cgMLSTv1.npz";
            my $precomputed_clusters_path_hcc = $precomputed_clusters_dir . lc($dir_name) . ".cgMLSTv1.HierCC.gz";
            my $local_ref_clusters_path_npz       = "$tmp_dir/precomputed_clusters.npz";
            my $local_ref_clusters_path_hcc       = "$tmp_dir/precomputed_clusters.HierCC.gz";
            my $heircc_out                = "$tmp_dir/cluster";
            run_cmd(["cp", $precomputed_clusters_path_npz, $local_ref_clusters_path_npz], "Copy Master npz");
            run_cmd(["cp", $precomputed_clusters_path_hcc, $local_ref_clusters_path_hcc], "Copy Master HierCC");

            # Load all existing cluster IDs from the reference HierCC file as a set
            my %ref_cluster_ids;
            if ($opt->dry_run) {
                print STDERR "DRY RUN - would load reference cluster IDs from $local_ref_clusters_path_hcc\n";
            } else {
                open(my $ref_fh, '-|', 'gzip', '-dc', $local_ref_clusters_path_hcc)
                    or die "Cannot open reference HierCC $local_ref_clusters_path_hcc: $!";
                my $ref_hdr = <$ref_fh>;  # skip header row
                while (my $ref_line = <$ref_fh>) {
                    chomp $ref_line;
                    my @ref_cols = split(/\t/, $ref_line);
                    shift @ref_cols;  # remove genome ID (first column)
                    for my $val (@ref_cols) {
                        $ref_cluster_ids{$val} = 1 if defined $val && $val ne '';
                    }
                }
                close($ref_fh);
                print STDERR "Loaded " . scalar(keys %ref_cluster_ids) . " unique cluster IDs from reference.\n";
            }

            run_cmd([
                "pHierCC",
                "--profile", "$tmp_dir/master_joined_clean.tsv",
                "--output",  $heircc_out,
                "--append",  $local_ref_clusters_path_npz
            ], "Cluster");
            
            my $new_clusters_path_hcc = "$tmp_dir/cluster.HierCC.gz";

            run_cmd(["gunzip", $new_clusters_path_hcc], "gunzip");

            my $heircc_unzipped = "$tmp_dir/cluster.HierCC";
            if ($opt->dry_run) {
                my %cgmlst_hc;
                $cgmlst_hc{"cgmlst_hc$_"} = "DRY_RUN" for (0, 2, 5, 10, 20, 50, 100);
                $sequence_typing->{cgmlst_hc} = \%cgmlst_hc;
                push(@{$gto_in->{sequence_types}}, $sequence_typing);
            } else {
                open(my $heircc_fh, '<', $heircc_unzipped) or die "Cannot open $heircc_unzipped: $!";
                my $heircc_header = <$heircc_fh>;
                chomp $heircc_header;
                my @heircc_keys = split(/\t/, $heircc_header);

                my ($heircc_data, $last_line);
                while (my $line = <$heircc_fh>) {
                    chomp $line;
                    my @cols = split(/\t/, $line);
                    $heircc_data = $line if $cols[0] eq 'input';
                    $last_line = $line;
                }
                close($heircc_fh);

                die "Could not find 'input' row in $heircc_unzipped\n" unless defined $heircc_data;

                if ($heircc_data ne $last_line) {
                    print STDERR "WARNING: 'input' row is not the last row in $heircc_unzipped\n";
                }

                my @heircc_values = split(/\t/, $heircc_data);
                shift @heircc_values;  # remove genome ID column
                shift @heircc_keys;    # remove genome ID header

                my %cgmlst_hc;
                for my $i (0 .. $#heircc_keys) {
                    my $k = $heircc_keys[$i];
                    my $v = $heircc_values[$i];

                    next unless $k =~ /^HC(\d+)$/i;
                    my $level = $1;
                    next unless $level =~ /^(0|2|5|10|20|50|100)$/;

                    if (exists $ref_cluster_ids{$v}) {
                        $cgmlst_hc{"cgmlst_hc$level"} = $v;
                    } else {
                        print STDERR "New cluster ID '$v' at HC$level not found in reference set - assigning null\n";
                        $cgmlst_hc{"cgmlst_hc$level"} = "null";
                    }
                }

                $sequence_typing->{cgmlst_hc} = \%cgmlst_hc;
                push(@{$gto_in->{sequence_types}}, $sequence_typing);
            }
        }

} else {
    print STDERR "Schema does not exist for this species\n";
}

if ($opt->out)
{
    $gto_in->destroy_to_file($opt->out);
}
else
{
    $gto_in->destroy_to_file(\*STDOUT);
}



sub generateShell{
    my ($output_shell, $content, $outshell ,$finish_string) = @_;
    unlink glob "$output_shell.*";
    $finish_string ||= "Still_waters_run_deep";
    chomp $content;
    open OUT,">$output_shell" or die "Cannot open file $output_shell:$!";
    print OUT "#!/bin/bash\n";
    print OUT "echo ==========start at : `date` ==========\n";
    print OUT "set -e \n";
    print OUT "$content && \\\n";
    print OUT "echo ==========end at : `date` ========== && \\\n";
    print OUT "echo $finish_string 1>&2 && \\\n";
    print OUT "echo $finish_string > $output_shell.sign\n";
    close OUT;
    my $run = "bash $output_shell 1>$output_shell.e 2>$output_shell.o";
    push @{$outshell}, $run;
}

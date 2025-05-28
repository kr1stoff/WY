#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Std;
use Cwd;

my($qsub_opt, %allJobs, $qsubDir, $shell, $wrn, $logf, $errorf, $finishedf);

use vars qw($opt_d $opt_b $opt_m $opt_s $opt_l $opt_r $opt_h);
getopts("d:b:m:s:l:rh");

if(defined $opt_h or @ARGV == 0){
    &usage();
    exit;
}

# 当前用户
chomp(my $user = `whoami`);

# 目录$qsubDir, 用于存放任务输出信息等
$wrn = "qsub_".time();
$qsubDir = $opt_d || $wrn;
my $opt_N = "$qsubDir/qsub";
$shell = shift;
$logf = "$qsubDir/log";
$errorf = "$qsubDir/error";
$finishedf = "$qsubDir/finished";
$qsub_opt = "/usr/bin/bash";
die "error opt_d!" if $qsubDir eq "/";
my $psu = `ps u -u $user`;
while ($psu =~ /$qsub_opt $opt_N\_(\d+)\.sh/gs) {
    print STDERR "error exists $psu";
    `sleep 10`;
}

# 默认每个sh文本放1个命令, 最大同时任务2
# 每隔60秒扫描任务状态, 最大尝试投递次数10
my $lines = $opt_b || 1;
my $maxJob = $opt_m || 2;
my $sleepTime = $opt_s || 60;
my $max_try = 10;
my $delay = $opt_l if defined $opt_l;
$max_try = 1 if (!$opt_r);
my $unfinished_number;

`mkdir -p $qsubDir`;
# 根据$shell文档生成并行运行的命令文档
my $split_number;
open IS, $shell or die $!;
while (<IS>) {
    chomp;
    $split_number++;
    my $num = 1;
    open OUTS, ">$opt_N\_$split_number.sh" or die $!;
    print OUTS "set -e\n";
    print OUTS $_;
    while ($num < $lines) {
        $num++;
        last if (eof(IS));
        chomp(my $command = <IS>);
        print OUTS "\n$command";
    }
    print OUTS "\necho this-work-is-complete\n";
    close OUTS;
    $allJobs{$split_number} = 1;
}
close IS;

&qsub_and_wait();

# 如果最大允许同时投递任务数目$maxJob小于总投递任务, 则投递$maxJob个任务
# 每隔$sleepTime时间查看任务状态, 如果在跑任务数目小于$sub_num, 则继续投递任务
# 将所有任务的运行结果写入$logf
sub qsub_and_wait {
    if (-e $finishedf) {
        open I, $finishedf;
        while (<I>) {
            if (/$opt_N\_(\d+)\.sh/gs) {
                delete $allJobs{$1};
            }
        }
        close I;
    }
    my @wait = sort {$a <=> $b} keys %allJobs;

    my $sub_num = $maxJob > @wait ? @wait : $maxJob;

    my $qnum = 0;
    my(%runJob, %error);
    while (@wait and $qnum < $sub_num) {
        my $i = shift @wait;
        sleep($delay) if defined $delay and $qnum > 0;
        my $w = "$qsub_opt $opt_N\_$i.sh > $opt_N\_$i.sh.o 2> $opt_N\_$i.sh.e &\n";
        print $w;
        `$w`;
        $runJob{$i} = "$opt_N\_$i.sh";
        $qnum++;
    }

    while (@wait or keys %runJob) {
        sleep($sleepTime);
        &check_job($user, \%error, \@wait, \%runJob);
        $qnum = keys %runJob;
        my $qnum0 = $qnum;
        while (@wait and $qnum < $sub_num) {
            my $i = shift @wait;
            sleep($delay) if defined $delay and $qnum > $qnum0;
            my $w = "$qsub_opt $opt_N\_$i.sh > $opt_N\_$i.sh.o 2> $opt_N\_$i.sh.e &\n";
            print $w;
            `$w`;
            $runJob{$i} = "$opt_N\_$i.sh";
            $qnum++;
        }
    }

    open OUTL, ">>$logf" or die $!;
    if (keys %error) {
        print OUTL "There are some job can't run finish, check the shell and qsub again\n";
        for (sort keys %error) {
            print OUTL "$_\n";
        }
        die "There are some job can't run finish, check the shell and qsub again";
    } else {
        print OUTL "All jobs are finished correctly\n";
    }
    close OUTL;
}

# 检查投递任务的状态, 运行ps u, 获得当前并行任务的ID
# 对于已经停止的任务, 如果是完成了, 则从正在跑的任务名单剔除, 不然且在错误次数少于最大限度时,重新加入等待名单
sub check_job {
    my($userName, $error, $wait, $run) = @_;
    my %running;
    my $psu = `ps u -u $userName`;
    while ($psu =~ /$qsub_opt $opt_N\_(\d+)\.sh/gs) {
        my $num = $1;
        next unless defined($$run{$num});
        $running{$num} = undef;
    }

    foreach my $id (sort {$a <=> $b} keys %$run) {
        my $split_shell = $$run{$id};
        if (!exists $running{$id}) {
            delete $$run{$id};
            chomp(my $log = `tail -1 $split_shell.o`);
            if ($log eq "this-work-is-complete") {
                delete($$error{$split_shell});
                `echo $split_shell is finished! >> $finishedf`;
            } else {
                `echo $split_shell has not finished! >> $errorf`;
                $$error{$split_shell}++;
                if ($$error{$split_shell} < $max_try) {
                    `rm $split_shell.[oe]`;
                    my $num = $1 if ($split_shell =~ /$opt_N\_(\d+)\.sh/);
                    unshift @$wait, $num;
                    `echo $split_shell has been reqsub >> $errorf`;
                }
            }
        }
    }
}

# 输出帮助信息
sub usage {
    print <<EOD;
usage: perl $0 [options] shell.sh
支持断点重跑, 需设置相同-d参数
    Options:
        -d  qsub script and log dir, default ./qsub_time()  # 工作目录, 存在则继续未完成任务, 保存拆分的小任务以及任务的输出信息等
        -b  set number of lines to form a job, default 1  # 指定每个子任务由多少行原脚本命令组成
        -m  set the maximum number of jobs to throw out, default 2  # 指定可同时提交任务的最大数目
        -s  set interval time of checking by qstat, default 60 seconds  # 指定每隔多长时间检查任务的状态
        -l  set delay time of each job, default no delay  # 指定投递任务间隔时间
        -r  mark to reqsub the job which was finished error, max reqsub 10 times, default not  # # 指定每个sh任务错误是否重投10次
        -h  show this help  # 显示帮助信息
EOD
}

package com.kmlab.util;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CondaEnvExecutor {
    private static final Logger logger = LogManager.getLogger(CondaEnvExecutor.class);
    private static final String condaActivate = getCondaActivatePath();

    public static String getCondaActivatePath() {
        Map<String, String> condaEnvProperties = PropertyLoader.loadCondaEnvProperties();
        String condaActivate = condaEnvProperties.get("conda.activate.path");
        return condaActivate;
    }

    // 方法重构, 添加默认参数, 默认不生成log文件
    public static void executeCommand(String command, String condaEnvName) {
        executeCommand(command, condaEnvName, null);
    }

    /**
     * 执行指定的命令，可以在指定的conda环境中运行，并将日志输出到文件。
     *
     * @param command      需要执行的命令字符串。
     * @param condaEnvName 指定的conda环境名称。
     * @param logFile      用于记录命令执行日志的文件路径，如果为空或null，则不记录日志。
     */
    public static void executeCommand(String command, String condaEnvName, String logFile) {
        List<String> safeCommand = new ArrayList<>();
        safeCommand.add("bash");
        safeCommand.add("-c");
        safeCommand.add(String.format("source %s %s && %s && conda deactivate", condaActivate, condaEnvName, command));
        logger.info("运行命令: " + String.join(" ", safeCommand));

        ProcessBuilder processBuilder = new ProcessBuilder(safeCommand);

        try {
            Process process = processBuilder.start();

            if (logFile != null && !logFile.isEmpty()) {
                try (PrintWriter stdLogWriter = new PrintWriter(new FileWriter(logFile))) {
                    Thread stdOutThread = new LogThread(process.getInputStream(), stdLogWriter, "stdout");
                    Thread stdErrThread = new LogThread(process.getErrorStream(), stdLogWriter, "stderr");
                    stdOutThread.start();
                    stdErrThread.start();

                    stdOutThread.join();
                    stdErrThread.join();
                } catch (IOException | InterruptedException e) {
                    e.printStackTrace();
                    throw new RuntimeException("日志写入或线程处理异常", e);
                }
            }

            int exitCode = process.waitFor();
            if (exitCode != 0) {
                throw new InterruptedException("命令执行失败，退出码: " + exitCode);
            }
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
            throw new RuntimeException("命令执行异常", e);
        }
    }

    private static class LogThread extends Thread {
        private final InputStream inputStream;
        private final PrintWriter writer;
        private final String logPrefix;

        public LogThread(InputStream inputStream, PrintWriter writer, String logPrefix) {
            this.inputStream = inputStream;
            this.writer = writer;
            this.logPrefix = logPrefix;
        }

        @Override
        public void run() {
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    writer.println(logPrefix + ": " + line);
                }
            } catch (IOException e) {
                e.printStackTrace();
                throw new RuntimeException("读取命令输出异常", e);
            }
        }
    }

    public static void parallelExecuteCommand(List<String[]> listCommandCondaLog, int parallelNumber) {
        ExecutorService executorService = Executors.newFixedThreadPool(parallelNumber);

        List<Future<Void>> futures = new ArrayList<>();
        for (String[] cmdDetails : listCommandCondaLog) {
            String command = cmdDetails[0];
            String condaEnvName = cmdDetails[1];
            String logFile = cmdDetails[2];

            Future<Void> future = executorService.submit(() -> {
                try {
                    executeCommand(command, condaEnvName, logFile);
                } catch (Exception e) {
                    e.printStackTrace();
                }
                return null;
            });
            futures.add(future);
        }
        executorService.shutdown();

        for (Future<Void> future : futures) {
            try {
                future.get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }

        while (!executorService.isTerminated()) {
            // 可能做一些清理工作或等待
        }

        logger.info("所有命令执行完成");
    }
}

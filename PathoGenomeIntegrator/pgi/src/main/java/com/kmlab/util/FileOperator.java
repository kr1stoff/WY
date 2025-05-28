package com.kmlab.util;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.regex.Pattern;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FileOperator {
    private static final Logger logger = LogManager.getLogger(FileOperator.class);

    /**
     * 在指定的位置创建指向源文件的软链接。
     * 如果目标路径已经存在软链接，则不会重新创建，只会记录日志。
     *
     * @param sourcePath     源文件的路径，将要指向这个路径的文件创建软链接。
     * @param linkTargetPath 软链接的目标路径，即软链接将要创建的位置。
     */
    public static void linkFile(Path sourcePath, Path linkTargetPath) {
        if (Files.exists(linkTargetPath)) {
            logger.debug("软链接已存在 " + linkTargetPath);
        } else {
            try {
                Files.createSymbolicLink(linkTargetPath, sourcePath);
            } catch (Exception e) {
                e.printStackTrace();
            }
            logger.debug(String.format("软链接创建成功 source: %s target: %s", linkTargetPath, sourcePath));
        }
    }

    /**
     * 在指定目录下查找包含特定内容且扩展名匹配的文件。
     *
     * @param directoryPath 指定的目录路径
     * @param targetContent 要在文件内容中查找的目标字符串
     * @param extensions    要查找的文件扩展名数组
     * @return 找到的第一个文件的路径，如果没有找到则返回空字符串
     */
    public static String findFileByContentAndExtensions(String directoryPath, String targetContent,
            String[] extensions) {
        // 将目录路径转换为Path对象
        Path startDir = Paths.get(directoryPath);

        // 构建正则表达式，用于匹配指定扩展名的文件
        StringBuilder extensionPattern = new StringBuilder();
        for (String extension : extensions) {
            if (extensionPattern.length() > 0) {
                extensionPattern.append("|");
            }
            extensionPattern.append(Pattern.quote(extension));
        }
        String regex = ".*" + Pattern.quote(targetContent) + ".*\\." + extensionPattern.toString() + "$";
        Pattern pattern = Pattern.compile(regex);
        // 使用AtomicReference来允许"修改"matchedFilePath
        AtomicReference<String> matchedFilePathRef = new AtomicReference<>("");

        try {
            // 使用FileVisitor遍历目录，查找匹配的文件
            SimpleFileVisitor<Path> visitor = new SimpleFileVisitor<Path>() {
                @Override
                public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                    // 如果文件名匹配正则表达式，则终止遍历
                    if (pattern.matcher(file.getFileName().toString()).matches()) {
                        matchedFilePathRef.set(file.toAbsolutePath().toString()); // 保存匹配文件的绝对路径
                        return FileVisitResult.TERMINATE;
                    }
                    return FileVisitResult.CONTINUE;
                }

                @Override
                public FileVisitResult visitFileFailed(Path file, IOException exc) {
                    // 忽略访问文件失败的情况
                    return FileVisitResult.CONTINUE;
                }
            };

            // 开始遍历目录
            Files.walkFileTree(startDir, visitor);
        } catch (IOException e) {
            // 处理遍历过程中的IO异常
            System.err.println("搜索过程中发生错误: " + e.getMessage());
        }
        return matchedFilePathRef.get(); // 返回找到的文件的绝对路径，如果没有找到则返回空字符串
    }

    /**
     * 递归方法用于查找具有指定扩展名的所有文件
     *
     * @param directory  要遍历的目录
     * @param extensions 文件后缀数组
     * @return 包含所有匹配文件的列表
     */
    public static List<File> findFilesWithExtensions(String directory, String[] extensions) {
        File fileDirectory = new File(directory);
        File[] fileList = fileDirectory.listFiles();
        List<File> filesWithExtensions = new ArrayList<>();

        if (fileList != null) {
            for (File file : fileList) {
                if (file.isDirectory()) {
                    // 如果是目录，则递归调用自身
                    filesWithExtensions.addAll(findFilesWithExtensions(file.toString(), extensions));
                } else {
                    // 检查文件后缀是否匹配
                    for (String ext : extensions) {
                        if (file.getName().toLowerCase().endsWith(ext.toLowerCase())) {
                            filesWithExtensions.add(file);
                            break; // 找到匹配后缀，跳出循环
                        }
                    }
                }
            }
        }
        return filesWithExtensions;
    }

    /**
     * 创建多个目录。
     *
     * @param parentDirectories 父目录的路径，必须为非空字符串。
     * @param sonDirectories    要创建的子目录数组，每个元素代表一个子目录。
     *                          说明：此方法会迭代遍历targetDirs数组，并为每个子目录在指定的parentDirectories下创建目录。
     *                          不会检查目录是否已存在，也不会处理任何异常情况，如权限不足等。
     */
    public static void mkdirs(String parentDirectories, String[] sonDirectories) {
        for (String sonDirectory : sonDirectories) {
            mkdir(parentDirectories + "/" + sonDirectory);
        }
    }

    /**
     * 创建目录。
     * 该方法会遍历传入的字符串数组，并对每个字符串代表的目录路径调用{@code mkdir}方法进行创建。
     *
     * @param targetDirs 目标目录数组，每个元素代表一个待创建的目录路径。
     *                   如果路径中包含父目录，也会一并创建。
     */
    public static void mkdirs(String[] targetDirs) {
        // 遍历目标目录数组，逐个创建目录
        for (String targetDir : targetDirs) {
            mkdir(targetDir); // 调用mkdir方法创建单个目录
        }
    }

    /**
     * 创建目录
     *
     * @param targetDir 指定要创建的目录路径
     *                  该方法尝试创建指定路径的目录。如果目录结构中的任一部分不存在，都将被创建。
     *                  如果出现IO异常，例如没有权限创建目录或磁盘空间不足，将打印异常堆栈跟踪。
     */
    public static void mkdir(String targetDir) {
        logger.debug("创建目录: " + targetDir);
        try {
            Path pathTargetDir = Paths.get(targetDir);
            Files.createDirectories(pathTargetDir);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}

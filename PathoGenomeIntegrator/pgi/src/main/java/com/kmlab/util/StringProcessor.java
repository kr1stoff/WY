package com.kmlab.util;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class StringProcessor {
    private static final Logger logger = LogManager.getLogger(StringProcessor.class);

    /**
     * 将输入字符串的第一个字母转换为大写。
     *
     * @param inputString 输入的字符串，要求非空。
     * @return 转换后的字符串，其中第一个字母将被转换为大写形式。
     */
    public static String toTitle(String inputString) {
        logger.info("将输入字符串的第一个字母转换为大写");
        char firstLetter = inputString.charAt(0);
        char upperFirstLetter = Character.toUpperCase(firstLetter);
        return inputString.replace(inputString.charAt(0), upperFirstLetter);
    }
}

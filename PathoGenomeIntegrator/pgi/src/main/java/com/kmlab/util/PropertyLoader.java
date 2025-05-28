package com.kmlab.util;

import java.util.Map;
import java.util.Properties;
import java.io.IOException;
import java.util.HashMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PropertyLoader {
    private static final Logger logger = LogManager.getLogger(PropertyLoader.class);
    private static final String BASIC_PROPERTY = "basic.properties";
    private static final String PARAMETERS_PROPERTY = "pararmeters.properties";
    private static final String CONDA_ENV_PROPERTY = "conda_env.properties";

    public static Map<String, String> propertyToMap(String propertyFilePath) {
        Properties properties = new Properties();
        Map<String, String> propertiesMap = new HashMap<>();

        try (var inputStream = PropertyLoader.class.getClassLoader().getResourceAsStream(propertyFilePath)) {
            if (inputStream != null) {
                properties.load(inputStream);
                for (String key : properties.stringPropertyNames()) {
                    propertiesMap.put(key, properties.getProperty(key));
                }
            } else {
                throw new RuntimeException("找不到配置文件: " + propertyFilePath);
            }
        } catch (IOException e) {
            throw new RuntimeException("加载配置文件失败: " + e.getMessage());
        }

        return propertiesMap;
    }

    public static Map<String, String> loadBasicProperties() {
        logger.info("加载基础配置");
        return propertyToMap(BASIC_PROPERTY);
    }

    public static Map<String, String> loadParametersProperties() {
        logger.info("加载参数配置");
        return propertyToMap(PARAMETERS_PROPERTY);
    }

    public static Map<String, String> loadCondaEnvProperties() {
        logger.info("加载Conda环境配置");
        return propertyToMap(CONDA_ENV_PROPERTY);
    }
}

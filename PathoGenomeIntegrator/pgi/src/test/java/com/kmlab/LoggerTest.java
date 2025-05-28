package com.kmlab;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;

public class LoggerTest {
    // private static final Log logger = LogFactory.getLog(LoggerTest.class);
    private static final Logger logger = LogManager.getLogger(LoggerTest.class);

    @Test
    public void test() {
        logger.debug("调试信息");
        logger.info("信息信息");
        logger.warn("警告信息");
        logger.error("错误信息");
    }
}

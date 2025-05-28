package com.kmlab.util;

import com.sun.management.OperatingSystemMXBean;
import java.lang.management.ManagementFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SystemInfo {
    private static final Logger logger = LogManager.getLogger(SystemInfo.class);

    public static int getThreadNumber() {
        OperatingSystemMXBean osBean = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();
        int availableProcessors = osBean.getAvailableProcessors();
        logger.info("系统最大可用线程数: " + availableProcessors);
        return availableProcessors;
    }
}

package com.kmlab.util;

import org.junit.Test;

public class SystemInfoTest {
    @Test
    public void testGetThreadNumber() {
        int availableProcessors = SystemInfo.getThreadNumber();
        System.out.println("availableProcessors: " + availableProcessors);
    }
}

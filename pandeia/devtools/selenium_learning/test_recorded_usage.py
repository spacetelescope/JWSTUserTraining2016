# -*- coding: utf-8 -*-
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import NoAlertPresentException
import unittest, time, re

class TestWorkbookPageElements(unittest.TestCase):
    def setUp(self):
        self.driver = webdriver.Firefox()
        self.driver.implicitly_wait(30)
        self.base_url = "http://etcbrady.stsci.edu:4242"
        self.verificationErrors = []
        self.accept_next_alert = True
    
    def test_workbook_page_elements(self):
        driver = self.driver
        driver.get(self.base_url + "/workbook.html")
        driver.find_element_by_link_text("Backgrounds").click()
        driver.find_element_by_link_text("Instrument Setup").click()
        driver.find_element_by_link_text("Strategy").click()
        driver.find_element_by_link_text("Detector Setup").click()
        driver.find_element_by_link_text("Warnings").click()
        driver.find_element_by_link_text("Focal Plane Rate").click()
        driver.find_element_by_link_text("Target").click()
        driver.find_element_by_link_text("Detector").click()
        driver.find_element_by_id("scene_library_menu").click()
        driver.find_element_by_link_text("Editor").click()
        driver.find_element_by_css_selector("div.modal-footer > button.btn-small.btn-primary").click()
        driver.get_cookies()
    
    def is_element_present(self, how, what):
        try: self.driver.find_element(by=how, value=what)
        except NoSuchElementException, e: return False
        return True
    
    def is_alert_present(self):
        try: self.driver.switch_to_alert()
        except NoAlertPresentException, e: return False
        return True
    
    def close_alert_and_get_its_text(self):
        try:
            alert = self.driver.switch_to_alert()
            alert_text = alert.text
            if self.accept_next_alert:
                alert.accept()
            else:
                alert.dismiss()
            return alert_text
        finally: self.accept_next_alert = True

    def is_cookie_present(self):
        try:
            cookies = self.driver.get_cookies()
            assert (len(cookies.keys()),0)
        except:
            raise NoSuchElementException

    
    def tearDown(self):
        self.driver.quit()
        self.assertEqual([], self.verificationErrors)

if __name__ == "__main__":
    unittest.main()

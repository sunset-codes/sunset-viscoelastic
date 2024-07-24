#!/bin/bash

# modify the below to account for how many characters we need to remove. fields10000_XX is to become fields10000
for x in fields*; do mv "$x" "${x%?}"; done


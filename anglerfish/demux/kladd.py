
attrs = ["name",
"sequence_token",
"index",
"has_index",
"len_index",
"len_before_index",
"len_after_index",
"has_umi",
"len_umi",
"len_umi_before_index",
"len_umi_after_index",]

for attr in attrs:
    e = eval(f"self.{attr}")
    print(f"self.{attr} = {e}")
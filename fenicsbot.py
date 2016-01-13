import twitter
from twitter.error import TwitterError
import traceback
from help_tweets import help_dict
from time import sleep, time
from parser import parse, excise


class FEniCSbot(object):
    def __init__(self, twitter_api_object,
                 sleep_time=10, break_loop=lambda *_: False):
        self.api = twitter_api_object
        self.sleep_time = sleep_time
        self._last_check_id = 1
        self.break_loop = break_loop

    def print_status(self, message):
        """
        Helper function for printing a short status to standard out.
        Made into a function and not just a print statement so it's 
        easy to change in the future, I guess?

        :param message: Message to print.
        """

        print message

    def tweet_image(self, img_fn, tweet):
        """
        Tweets an image as a reply tweet.

        :param img_fn: Filename of image to tweet.
        :param tweet: Twitter status object of tweet to reply to.
        """

        self.print_status("Tweeting a solution!")
        solution_tweet_text = ".@{}: I solved your problem ({})!".format(tweet.user.screen_name, 
                                                                        excise(tweet.text).strip())[:100]
        self.api.PostMedia(solution_tweet_text, img_fn,
                           in_reply_to_status_id=str(tweet.id))

    def tweet_error(self, tweet, e):
        """
        Tweets an error message in reply to a tweet which did not parse.

        :param tweet: Tweet to reply to.
        :param e: The Exception object
        """
        self.print_status("Got error {} when replying to {}".format(e, str(tweet.text)))
        self.print_status("Traceback: ")
        traceback.print_exc()

        doc_desc = " See docs: http://bit.ly/1T4KuNt"

        error_tweet = ".@{}: I failed to solve your problem... {} {}"
        error_tweet = error_tweet.format(tweet.user.screen_name, 
                                         e.message, doc_desc)

        self.print_status(error_tweet)
        self.api.PostUpdate(error_tweet[:140], 
                            in_reply_to_status_id=str(tweet.id))

    def tweet_welcome(self, tweet):
        """
        Replies to a thank you tweet with something nice.
        
        :param tweet: Tweet to reply to.
        """
        welcome_message = ".@{}: You're welcome!".format(tweet.user.screen_name)
        self.api.PostUpdate(welcome_message[:140], 
                            in_reply_to_status_id=str(tweet.id))
    
    def tweet_help(self, tweet):
        """
        Replies to a tweet requesting help about a specific model/command.
        
        :param tweet: Tweet to reply to.
        """
        
        help_pos = tweet.text.find("help")
        help_subject = tweet.text[help_pos+4:].split()[0].lower()

        if help_subject in help_dict:
            help_message = help_dict[help_subject]
        else:
            help_message = help_dict["documentation"]
        
        self.api.PostUpdate(help_message[:140], 
                            in_reply_to_status_id=str(tweet.id))
        
    
    
    def get_mentions(self):
        """
        Gets tweets @FEniCSbot, and returns those which are recent, have not
        been answered yet, are not from @FEniCSbot, and include the word 'solve'.
        
        Returns these tweets categorized into solve requests, thank you messages
        and requests for help.
        """

        new_mentions = self.api.GetSearch(term="@fenicsbot",
                                          since_id=self._last_check_id)
        if self._last_check_id == 1:
                # in first iteration, don't grab tweets from beginning of time
                is_recent = lambda t: (t.created_at_in_seconds >
                                       self.program_start_time)
                new_mentions = filter(lambda t: is_recent(t), new_mentions)

        if len(new_mentions) > 0:
            self._last_check_id = new_mentions[0].id


        # list of predicates which are True for solve requests
        solve_criteria = [
            # ignore tweets from self
            lambda t: t.user.screen_name.lower() != "fenicsbot",

            # FEniCSbot should not reply to tweets not including "solve"
            lambda t: "solve" in t.text.lower(),

            # # FEniCSbot should not reply to retweets
            # lambda t: not hasattr(t, "retweeted_status")
            lambda t: t.text[:3] != "RT "
        ]
        
        thanks_criteria = [lambda t: "thank" in t.text.lower()]
        help_criteria = [lambda t: "help" in t.text.lower()]

        list_of_category_criterion_lists = [
            solve_criteria, thanks_criteria, help_criteria
        ]
    
        categories = []
        for category_criterion_list in list_of_category_criterion_lists:
            categories.append(new_mentions)
            for criterion in category_criterion_list:
                categories[-1] = filter(criterion, categories[-1])
            new_mentions = filter(lambda t: t not in categories[-1], new_mentions)
        return categories


    def start(self):
        """
        Starts FEniCSbot listen loop, causing it to
        reply to all tweets containing @FEniCSbot.
        """
        self.program_start_time = time()

        while not self.break_loop(self):
            self.print_status("\nScanning for new tweets...")
            
            solve_requests, thank_yous, help_requests = self.get_mentions()
            for tweet in solve_requests:
                try:
                    solver = parse(tweet.text)
                    solver.solve()
                    img_fn = solver.plot()
                    self.tweet_image(img_fn, tweet)

                except Exception as e :
                    self.tweet_error(tweet, e)
                        
            for tweet in thank_yous:
                self.tweet_welcome(tweet)
            
            for tweet in help_requests:
                self.tweet_help(tweet)

            self.print_status("Scan for new tweets complete.\n")
            sleep(self.sleep_time)




if __name__ == "__main__":
    from secrets import secret_dict

    twitter_api = twitter.Api(**secret_dict)
    
    while True:
        try:
            bot = FEniCSbot(twitter_api)
            bot.start()
        except TwitterError as te:
            continue
